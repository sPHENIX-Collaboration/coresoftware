#include "RetowerCEMC.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// standard includes
#include <cstdlib>
#include <iostream>
#include <utility>

RetowerCEMC::RetowerCEMC(const std::string &name)
  : SubsysReco(name)
{
}

int RetowerCEMC::InitRun(PHCompositeNode *topNode)
{
  CreateNode(topNode);
  get_first_phi_index(topNode);
  if (_weighted_energy_distribution == 1)
  {
    get_weighted_fraction(topNode);
  }
  else
  {
    get_fraction(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "RetowerCEMC::process_event: entering" << std::endl;
  }

  RawTowerContainer *towersEM3 = nullptr;
  TowerInfoContainer *towerinfosEM3 = nullptr;
  if (m_use_towerinfo)
  {
    EMTowerName = m_towerNodePrefix + "_CEMC";
    towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
    if (!towerinfosEM3)
    {
      std::cout << "RetowerCEMC::process_event: Cannot find node " << EMTowerName << std::endl;
      exit(1);
    }
  }
  else
  {
    towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;
    }
  }

  if (m_use_towerinfo)
  {
    unsigned int nchannels = towerinfosEM3->size();
    for (unsigned int channel = 0; channel < nchannels; channel++)
    {
      TowerInfo *tower = towerinfosEM3->get_tower_at_channel(channel);
      unsigned int channelkey = towerinfosEM3->encode_key(channel);
      int ieta = towerinfosEM3->getTowerEtaBin(channelkey);
      int iphi = towerinfosEM3->getTowerPhiBin(channelkey);
      rawtower_e[ieta][iphi] = tower->get_energy();
      rawtower_time[ieta][iphi] = tower->get_time_float();
      rawtower_status[ieta][iphi] = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2();
    }
    EMRetowerName = m_towerNodePrefix + "_CEMC_RETOWER";
    TowerInfoContainer *emcal_retower = findNode::getClass<TowerInfoContainer>(topNode, EMRetowerName);
    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: filling " << EMRetowerName << " node" << std::endl;
    }
    for (int ieta_ihcal = 0; ieta_ihcal < neta_ihcal; ++ieta_ihcal)
    {
      for (int iphi_ihcal = 0; iphi_ihcal < nphi_ihcal; ++iphi_ihcal)
      {
        double retower_e_temp = 0;
        double retower_time_temp = 0;
        double retower_badarea = 0;
        for (int ieta_emcal = retower_lowerbound_originaltower_ieta[ieta_ihcal]; ieta_emcal <= retower_upperbound_originaltower_ieta[ieta_ihcal]; ++ieta_emcal)
        {
          for (int iphi_emcal = retower_first_lowerbound_originaltower_iphi + (iphi_ihcal * 4); iphi_emcal < retower_first_lowerbound_originaltower_iphi + iphi_ihcal * 4 + 4; ++iphi_emcal)
          {
            int iphi_emcal_wrap = iphi_emcal;
            if (iphi_emcal > nphi_emcal - 1)
            {
              iphi_emcal_wrap -= nphi_emcal;
            }
            double fraction_temp;
            if (ieta_emcal == retower_lowerbound_originaltower_ieta[ieta_ihcal])
            {
              fraction_temp = retower_lowerbound_originaltower_fraction[ieta_ihcal];
            }
            else if (ieta_emcal == retower_upperbound_originaltower_ieta[ieta_ihcal])
            {
              fraction_temp = retower_upperbound_originaltower_fraction[ieta_ihcal];
            }
            else
            {
              fraction_temp = 1;
            }
            if (rawtower_status[ieta_emcal][iphi_emcal_wrap])
            {
              retower_badarea += fraction_temp;
            }
            else
            {
              retower_e_temp += rawtower_e[ieta_emcal][iphi_emcal_wrap] * fraction_temp;
              retower_time_temp += rawtower_time[ieta_emcal][iphi_emcal_wrap] * rawtower_e[ieta_emcal][iphi_emcal_wrap] * fraction_temp;
            }
          }
        }
        unsigned int towerkey = TowerInfoDefs::encode_hcal(ieta_ihcal, iphi_ihcal);
        unsigned int towerindex = emcal_retower->decode_key(towerkey);
        TowerInfo *towerinfo = emcal_retower->get_tower_at_channel(towerindex);
        double scalefactor = retower_badarea / retower_totalarea[ieta_ihcal];
        if (scalefactor > _frac_cut)
        {
          towerinfo->set_energy(0);
          towerinfo->set_isHot(true);
        }
        else
        {
          towerinfo->set_energy(retower_e_temp / (double) (1 - scalefactor));
          if (retower_e_temp == 0)
          {
            towerinfo->set_time_float(0);
          }
          else
          {
            towerinfo->set_time_float((retower_time_temp / retower_e_temp));
          }
          towerinfo->set_chi2(scalefactor);
        }
      }
    }
  }
  else
  {
    for (int ieta_ihcal = 0; ieta_ihcal < neta_ihcal; ++ieta_ihcal)
    {
      for (int iphi_ihcal = 0; iphi_ihcal < nphi_ihcal; ++iphi_ihcal)
      {
        double retower_e_temp = 0;
        for (int ieta_emcal = retower_lowerbound_originaltower_ieta[ieta_ihcal]; ieta_emcal <= retower_upperbound_originaltower_ieta[ieta_ihcal]; ++ieta_emcal)
        {
          for (int iphi_emcal = retower_first_lowerbound_originaltower_iphi; iphi_emcal < retower_first_lowerbound_originaltower_iphi + iphi_ihcal * 4; ++iphi_emcal)
          {
            int iphi_emcal_wrap = iphi_emcal;
            if (iphi_emcal > nphi_emcal - 1)
            {
              iphi_emcal_wrap -= nphi_emcal;
            }
            RawTower *tower = towersEM3->getTower(ieta_emcal, iphi_emcal_wrap);
            double energy = tower->get_energy();
            double fraction_temp;
            if (ieta_emcal == retower_lowerbound_originaltower_ieta[ieta_ihcal])
            {
              fraction_temp = retower_lowerbound_originaltower_fraction[ieta_ihcal];
            }
            else if (ieta_emcal == retower_upperbound_originaltower_ieta[ieta_ihcal])
            {
              fraction_temp = retower_upperbound_originaltower_fraction[ieta_ihcal];
            }
            else
            {
              fraction_temp = 1;
            }
            retower_e_temp += energy * fraction_temp;
          }
        }
        RawTowerContainer *emcal_retower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
        if (Verbosity() > 0)
        {
          std::cout << "RetowerCEMC::process_event: filling TOWER_CALIB_CEMC_RETOWER node, with initial size = " << emcal_retower->size() << std::endl;
        }
        RawTower *new_tower = new RawTowerv1();
        new_tower->set_energy(retower_e_temp);
        emcal_retower->AddTower(ieta_ihcal, iphi_ihcal, new_tower);
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << "RetowerCEMC::process_event: exiting" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHCompositeNode *emcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));
  if (!emcalNode)
  {
    std::cout << PHWHERE << "EMCal Node note found, doing nothing." << std::endl;
  }

  if (m_use_towerinfo)
  {
    EMRetowerName = m_towerNodePrefix + "_CEMC_RETOWER";
    IHTowerName = m_towerNodePrefix + "_HCALIN";
    TowerInfoContainer *test_emcal_retower = findNode::getClass<TowerInfoContainer>(topNode, EMRetowerName);
    TowerInfoContainer *hcal_towers = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
    if (!test_emcal_retower)
    {
      if (Verbosity() > 0)
      {
        std::cout << "RetowerCEMC::CreateNode : creating " << EMRetowerName << " node " << std::endl;
      }
      if (!hcal_towers)
      {
        std::cout << PHWHERE << " Could not find input HCAL tower node: " << IHTowerName << std::endl;
        exit(1);
      }
      TowerInfoContainer *emcal_retower = dynamic_cast<TowerInfoContainer *>(hcal_towers->CloneMe());
      PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(emcal_retower, EMRetowerName, "PHObject");
      emcalNode->addNode(emcalTowerNode);
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "RetowerCEMC::CreateNode : " << EMRetowerName << " already exists! " << std::endl;
      }
    }
  }
  else
  {
    RawTowerContainer *test_emcal_retower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
    if (!test_emcal_retower)
    {
      if (Verbosity() > 0)
      {
        std::cout << "RetowerCEMC::CreateNode : creating TOWER_CALIB_CEMC_RETOWER node " << std::endl;
      }
      RawTowerContainer *emcal_retower = new RawTowerContainer(RawTowerDefs::CalorimeterId::HCALIN);
      PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(emcal_retower, "TOWER_CALIB_CEMC_RETOWER", "PHObject");
      emcalNode->addNode(emcalTowerNode);
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "RetowerCEMC::CreateNode : TOWER_CALIB_CEMC_RETOWER already exists! " << std::endl;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void RetowerCEMC::get_first_phi_index(PHCompositeNode *topNode)
{
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!geomEM || !geomIH)
  {
    std::cout << "RetowerCEMC::get_first_phi_index: Could not locate geometry nodes" << std::endl;
    exit(1);
  }

  bool find_first_lowerbound = false;
  int iphi_emcal = 0;
  while (iphi_emcal < nphi_emcal)
  {
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, iphi_emcal);
    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(key);
    int this_IHphibin = geomIH->get_phibin(tower_geom->get_phi());
    if (this_IHphibin == 0)
    {
      find_first_lowerbound = true;
      break;
    }
    iphi_emcal++;
  }

  if (find_first_lowerbound && iphi_emcal == 0)
  {
    bool outofrange = false;
    int iphi_emcal_temp = nphi_emcal - 1;
    while (iphi_emcal_temp > iphi_emcal)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, iphi_emcal_temp);
      RawTowerGeom *tower_geom = geomEM->get_tower_geometry(key);
      int this_IHphibin = geomIH->get_phibin(tower_geom->get_phi());
      if (this_IHphibin == nphi_ihcal - 1)
      {
        outofrange = true;
        break;
      }
      iphi_emcal_temp--;
    }
    if (!outofrange)
    {
      std::cout << "RetowerCEMC::get_first_phi_index: could not find matching EMCal towers for iphi_ihcal = " << nphi_ihcal - 1 << std::endl;
      exit(1);
    }
    if (iphi_emcal_temp + 1 == nphi_emcal)
    {
      retower_first_lowerbound_originaltower_iphi = 0;
    }
    else
    {
      retower_first_lowerbound_originaltower_iphi = iphi_emcal_temp + 1;
    }
  }
  else if (!find_first_lowerbound)
  {
    std::cout << "RetowerCEMC::get_first_phi_index: could not find matching EMCal towers for iphi_ihcal = 0" << std::endl;
    exit(1);
  }
  else
  {
    retower_first_lowerbound_originaltower_iphi = iphi_emcal;
  }
}

void RetowerCEMC::get_fraction(PHCompositeNode *topNode)
{
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!geomEM || !geomIH)
  {
    std::cout << "RetowerCEMC::get_fraction: Could not locate geometry nodes" << std::endl;
    exit(1);
  }

  int ieta_emcal = 0;
  int ieta_ihcal = 0;
  while (ieta_emcal < neta_emcal && ieta_ihcal < neta_ihcal)
  {
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta_emcal, 0);
    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(key);
    int this_IHetabin = geomIH->get_etabin(tower_geom->get_eta());
    if (this_IHetabin == ieta_ihcal)
    {
      retower_lowerbound_originaltower_ieta[ieta_ihcal] = ieta_emcal;
      ieta_ihcal++;
    }
    ieta_emcal++;
  }
  for (int ieta = 0; ieta < neta_ihcal; ++ieta)
  {
    if (ieta < neta_ihcal - 1)
    {
      retower_upperbound_originaltower_ieta[ieta] = retower_lowerbound_originaltower_ieta[ieta + 1] - 1;
    }
    else
    {
      retower_upperbound_originaltower_ieta[ieta] = neta_emcal - 1;
    }
    retower_lowerbound_originaltower_fraction[ieta] = 1;
    retower_upperbound_originaltower_fraction[ieta] = 1;
  }
}

void RetowerCEMC::get_weighted_fraction(PHCompositeNode *topNode)
{
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!geomEM || !geomIH)
  {
    std::cout << "RetowerCEMC::get_weighted_fraction: Could not locate geometry nodes" << std::endl;
    exit(1);
  }

  int ieta_emcal = 0;
  for (int ieta_ihcal = 0; ieta_ihcal < neta_ihcal; ++ieta_ihcal)
  {
    std::pair<double, double> range_ihcal = geomIH->get_etabounds(ieta_ihcal);
    double ihcal_lowerbound = range_ihcal.first;
    double ihcal_upperbound = range_ihcal.second;

    bool found_lowerbound = false;
    bool found_upperbound = false;
    while ((!found_lowerbound || !found_upperbound) && ieta_emcal < neta_emcal)
    {
      std::pair<double, double> range_emcal = geomEM->get_etabounds(ieta_emcal);
      double emcal_lowerbound = range_emcal.first;
      double emcal_upperbound = range_emcal.second;
      if (!found_lowerbound)
      {
        if (emcal_upperbound > ihcal_lowerbound && emcal_lowerbound <= ihcal_lowerbound)
        {
          retower_lowerbound_originaltower_ieta[ieta_ihcal] = ieta_emcal;
          retower_lowerbound_originaltower_fraction[ieta_ihcal] = (emcal_upperbound - ihcal_lowerbound) / (emcal_upperbound - emcal_lowerbound);
          found_lowerbound = true;
        }
        if (emcal_upperbound > ihcal_lowerbound && emcal_lowerbound > ihcal_lowerbound)
        {
          retower_lowerbound_originaltower_ieta[ieta_ihcal] = ieta_emcal;
          retower_lowerbound_originaltower_fraction[ieta_ihcal] = 1;
          found_lowerbound = true;
        }
      }
      else
      {
        if (emcal_upperbound >= ihcal_upperbound && emcal_lowerbound < ihcal_upperbound)
        {
          retower_upperbound_originaltower_ieta[ieta_ihcal] = ieta_emcal;
          retower_upperbound_originaltower_fraction[ieta_ihcal] = (ihcal_upperbound - emcal_lowerbound) / (emcal_upperbound - emcal_lowerbound);
          found_upperbound = true;
        }
        if (emcal_upperbound > ihcal_upperbound && emcal_lowerbound > ihcal_upperbound)
        {
          ieta_emcal--;
          retower_upperbound_originaltower_ieta[ieta_ihcal] = ieta_emcal;
          retower_upperbound_originaltower_fraction[ieta_ihcal] = 1;
          found_upperbound = true;
        }
      }
      if (found_lowerbound && found_upperbound)
      {
        retower_totalarea[ieta_ihcal] = 4 * (retower_lowerbound_originaltower_fraction[ieta_ihcal] + retower_upperbound_originaltower_fraction[ieta_ihcal] + (retower_upperbound_originaltower_ieta[ieta_ihcal] - retower_lowerbound_originaltower_ieta[ieta_ihcal] - 1));
      }
      else
      {
        ieta_emcal++;
      }
    }
    if (!found_lowerbound)
    {
      std::cout << "RetowerCEMC::get_weighted_fraction: could not find lower bound EMCal towers for ieta_ihcal = " << ieta_ihcal << std::endl;
      exit(1);
    }
    else if (!found_upperbound)
    {
      std::cout << "RetowerCEMC::get_weighted_fraction: could not find upper bound EMCal towers for ieta_ihcal = " << ieta_ihcal << std::endl;
      exit(1);
    }
  }
}

#include "RetowerCEMC.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

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
#include <iostream>
#include <map>
#include <utility>
#include <vector>

RetowerCEMC::RetowerCEMC(const std::string &name)
  : SubsysReco(name)
{
}

int RetowerCEMC::InitRun(PHCompositeNode *topNode)
{
  CreateNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "RetowerCEMC::process_event: entering" << std::endl;
  }

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = nullptr;
  TowerInfoContainer *towerinfosEM3 = nullptr;
  if (m_use_towerinfo)
  {
    towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  }
  else
  {
    towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;
    }
  }

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  // setup grid

  if (_NETA < 0)
  {
    _NETA = geomIH->get_etabins();
    _NPHI = geomIH->get_phibins();

    _EMCAL_RETOWER_E.resize(_NETA, std::vector<float>(_NPHI, 0));
    _EMCAL_RETOWER_T.resize(_NETA, std::vector<int>(_NPHI, 0));
  }

  for (int eta = 0; eta < _NETA; eta++)
  {
    for (int phi = 0; phi < _NPHI; phi++)
    {
      _EMCAL_RETOWER_E[eta][phi] = 0;
      _EMCAL_RETOWER_T[eta][phi] = 0;
    }
  }

  // partition existing CEMC energies among grid
  if (m_use_towerinfo)
  {
    if (!towerinfosEM3)
    {
      std::cout << PHWHERE << "no emcal tower info object, doing nothing" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    unsigned int nchannels = towerinfosEM3->size();
    for (unsigned int channel = 0; channel < nchannels; channel++)
    {
      TowerInfo *tower = towerinfosEM3->get_tower_at_channel(channel);
      unsigned int channelkey = towerinfosEM3->encode_key(channel);
      int ieta = towerinfosEM3->getTowerEtaBin(channelkey);
      int iphi = towerinfosEM3->getTowerPhiBin(channelkey);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
      RawTowerGeom *tower_geom = geomEM->get_tower_geometry(key);
      int this_IHetabin = geomIH->get_etabin(tower_geom->get_eta());
      double fractionalcontribution[3] = {0};

      if (_WEIGHTED_ENERGY_DISTRIBUTION == 1)
      {
        std::pair<double, double> range_embin = geomEM->get_etabounds(tower_geom->get_bineta());
        for (int etabin_iter = -1; etabin_iter <= 1; etabin_iter++)
        {
          if (this_IHetabin + etabin_iter < 0 || this_IHetabin + etabin_iter >= _NETA)
          {
            continue;
          }
          std::pair<double, double> range_ihbin = geomIH->get_etabounds(this_IHetabin + etabin_iter);
          if (range_ihbin.first <= range_embin.first && range_ihbin.second >= range_embin.second)
          {
            fractionalcontribution[etabin_iter + 1] = 1;
          }
          else if (range_ihbin.first <= range_embin.first && range_ihbin.second < range_embin.second && range_embin.first < range_ihbin.second)
          {
            fractionalcontribution[etabin_iter + 1] = (range_ihbin.second - range_embin.first) / (range_embin.second - range_embin.first);
          }
          else if (range_ihbin.first > range_embin.first && range_ihbin.second >= range_embin.second && range_embin.second > range_ihbin.first)
          {
            fractionalcontribution[etabin_iter + 1] = (range_embin.second - range_ihbin.first) / (range_embin.second - range_embin.first);
          }
          else
          {
            fractionalcontribution[etabin_iter + 1] = 0;
          }
        }
      }
      else
      {
        fractionalcontribution[0] = 0;
        fractionalcontribution[1] = 1;
        fractionalcontribution[2] = 0;
      }

      int this_IHphibin = geomIH->get_phibin(tower_geom->get_phi());
      float this_E = tower->get_energy();
      int this_T = tower->get_time();
      for (int etabin_iter = -1; etabin_iter <= 1; etabin_iter++)
      {
        if (this_IHetabin + etabin_iter < 0 || this_IHetabin + etabin_iter >= _NETA)
        {
          continue;
        }
        _EMCAL_RETOWER_E[this_IHetabin + etabin_iter][this_IHphibin] += this_E * fractionalcontribution[etabin_iter + 1];
        if (this_T == -10 || this_T == -11)  // if a tower in this retower is masked, mask the retower as well. Masking is noted by time = -10 or -11
        {
          // currently just set the retowered tower to -10 even if the original tower was -11
          _EMCAL_RETOWER_T[this_IHetabin + etabin_iter][this_IHphibin] = -10;
        }
      }
    }
  }

  else
  {
    RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());

      int this_IHetabin = geomIH->get_etabin(tower_geom->get_eta());
      double fractionalcontribution[3] = {0};

      // distribute energy based on shadowing of the inner hcal geometry
      if (_WEIGHTED_ENERGY_DISTRIBUTION == 1)
      {
        std::pair<double, double> range_embin = geomEM->get_etabounds(tower_geom->get_bineta());
        for (int etabin_iter = -1; etabin_iter <= 1; etabin_iter++)
        {
          if (this_IHetabin + etabin_iter < 0 || this_IHetabin + etabin_iter >= _NETA)
          {
            continue;
          }
          std::pair<double, double> range_ihbin = geomIH->get_etabounds(this_IHetabin + etabin_iter);
          if (range_ihbin.first <= range_embin.first && range_ihbin.second >= range_embin.second)
          {
            fractionalcontribution[etabin_iter + 1] = 1;
          }
          else if (range_ihbin.first <= range_embin.first && range_ihbin.second < range_embin.second && range_embin.first < range_ihbin.second)
          {
            fractionalcontribution[etabin_iter + 1] = (range_ihbin.second - range_embin.first) / (range_embin.second - range_embin.first);
          }
          else if (range_ihbin.first > range_embin.first && range_ihbin.second >= range_embin.second && range_embin.second > range_ihbin.first)
          {
            fractionalcontribution[etabin_iter + 1] = (range_embin.second - range_ihbin.first) / (range_embin.second - range_embin.first);
          }
          else
          {
            fractionalcontribution[etabin_iter + 1] = 0;
          }
        }
      }
      else
      {
        fractionalcontribution[0] = 0;
        fractionalcontribution[1] = 1;
        fractionalcontribution[2] = 0;
      }

      int this_IHphibin = geomIH->get_phibin(tower_geom->get_phi());
      float this_E = tower->get_energy();

      for (int etabin_iter = -1; etabin_iter <= 1; etabin_iter++)
      {
        if (this_IHetabin + etabin_iter < 0 || this_IHetabin + etabin_iter >= _NETA)
        {
          continue;
        }
        _EMCAL_RETOWER_E[this_IHetabin + etabin_iter][this_IHphibin] += this_E * fractionalcontribution[etabin_iter + 1];
      }
    }
  }

  if (m_use_towerinfo)
  {
    TowerInfoContainer *emcal_retower = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: filling TOWERINFO_CALIB_CEMC_RETOWER node" << std::endl;
    }
    // create new towers
    for (int eta = 0; eta < _NETA; eta++)
    {
      for (int phi = 0; phi < _NPHI; phi++)
      {
        unsigned int towerkey = (((unsigned int) eta) << 16U) + phi;
        unsigned int towerindex = emcal_retower->decode_key(towerkey);
        TowerInfo *towerinfo = emcal_retower->get_tower_at_channel(towerindex);
        towerinfo->set_energy(_EMCAL_RETOWER_E[eta][phi]);
        if (_EMCAL_RETOWER_T[eta][phi] == -10)
        {
          towerinfo->set_energy(0);
        }
        towerinfo->set_time(_EMCAL_RETOWER_T[eta][phi]);
      }
    }
  }
  else
  {
    RawTowerContainer *emcal_retower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: filling TOWER_CALIB_CEMC_RETOWER node, with initial size = " << emcal_retower->size() << std::endl;
    }
    // create new towers
    for (int eta = 0; eta < _NETA; eta++)
    {
      for (int phi = 0; phi < _NPHI; phi++)
      {
        RawTower *new_tower = new RawTowerv1();

        new_tower->set_energy(_EMCAL_RETOWER_E[eta][phi]);
        emcal_retower->AddTower(eta, phi, new_tower);
      }
    }

    if (Verbosity() > 0)
    {
      std::cout << "RetowerCEMC::process_event: finished filling TOWER_CALIB_CEMC_RETOWER node, with final size = " << emcal_retower->size() << std::endl;
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

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the new EMCal towers
  PHCompositeNode *emcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));
  if (!emcalNode)
  {
    std::cout << PHWHERE << "EMCal Node note found, doing nothing." << std::endl;
  }

  if (m_use_towerinfo)
  {
    TowerInfoContainer *test_emcal_retower = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
    TowerInfoContainer *hcal_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (!test_emcal_retower)
    {
      if (Verbosity() > 0)
      {
        std::cout << "RetowerCEMC::CreateNode : creating TOWERINFO_CALIB_CEMC_RETOWER node " << std::endl;
      }

      TowerInfoContainer *emcal_retower = dynamic_cast<TowerInfoContainer *>(hcal_towers->CloneMe());
      PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(emcal_retower, "TOWERINFO_CALIB_CEMC_RETOWER", "PHObject");
      emcalNode->addNode(emcalTowerNode);
    }
    else
    {
      std::cout << "RetowerCEMC::CreateNode : TOWERINFO_CALIB_CEMC_RETOWER already exists! " << std::endl;
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
      std::cout << "RetowerCEMC::CreateNode : TOWER_CALIB_CEMC_RETOWER already exists! " << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

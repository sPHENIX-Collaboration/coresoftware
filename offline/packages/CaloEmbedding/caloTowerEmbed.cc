#include "caloTowerEmbed.h"

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
caloTowerEmbed::caloTowerEmbed(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "caloTowerEmbed::caloTowerEmbed(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
caloTowerEmbed::~caloTowerEmbed()
{
  std::cout << "caloTowerEmbed::~caloTowerEmbed() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int caloTowerEmbed::InitRun(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  PHCompositeNode *simTopNode = se->topNode("TOPSim");

  EventHeader *evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  if (evtHeader)
  {
    m_runNumber = evtHeader->get_RunNumber();
  }
  else
  {
    m_runNumber = 0;
  }
  if (Verbosity())
  {
    std::cout << "at run" << m_runNumber << std::endl;
  }

  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  topNode->print();
  simTopNode->print();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTowerEmbed::process_event(PHCompositeNode * /*topNode*/)
{
  ++m_eventNumber;

  for (int caloIndex = 0; caloIndex < 3; caloIndex++)
  {
    if (Verbosity())
    {
      std::cout << "event " << m_eventNumber << " working on caloIndex " << caloIndex;
    }
    if (caloIndex == 0)
    {
      std::cout << " CEMC" << std::endl;
    }
    else if (caloIndex == 1)
    {
      std::cout << " HCALIN" << std::endl;
    }
    else if (caloIndex == 2)
    {
      std::cout << " HCALOUT" << std::endl;
    }
    RawTowerDefs::keytype keyData = 0;
    RawTowerDefs::keytype keySim = 0;

    unsigned int ntowers = _data_towers[caloIndex]->size();
    for (unsigned int channel = 0; channel < ntowers; channel++)
    {
      unsigned int data_key = _data_towers[caloIndex]->encode_key(channel);
      unsigned int sim_key = _sim_towers[caloIndex]->encode_key(channel);

      int ieta_data = _data_towers[caloIndex]->getTowerEtaBin(data_key);
      int iphi_data = _data_towers[caloIndex]->getTowerPhiBin(data_key);
      int ieta_sim = _sim_towers[caloIndex]->getTowerEtaBin(sim_key);
      int iphi_sim = _sim_towers[caloIndex]->getTowerPhiBin(sim_key);

      if (caloIndex == 0)
      {
        if (!m_useRetower)
        {
          keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta_data, iphi_data);
          keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta_sim, iphi_sim);
        }
        else
        {
          keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_data, iphi_data);
          keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_sim, iphi_sim);
        }
      }
      else if (caloIndex == 1)
      {
        keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_data, iphi_data);
        keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_sim, iphi_sim);
      }
      else if (caloIndex == 2)
      {
        keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta_data, iphi_data);
        keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta_sim, iphi_sim);
      }

      TowerInfo *caloinfo_data = _data_towers[caloIndex]->get_tower_at_channel(channel);
      TowerInfo *caloinfo_sim = _sim_towers[caloIndex]->get_tower_at_channel(channel);

      float data_E = caloinfo_data->get_energy();
      float sim_E = caloinfo_sim->get_energy();
      float embed_E = data_E + sim_E;

      float data_phi = 0.0;
      float sim_phi = 0.0;

      float data_eta = 0.0;
      float sim_eta = 0.0;

      if (caloIndex == 0)
      {
        if (!m_useRetower)
        {
          data_phi = tower_geom->get_tower_geometry(keyData)->get_phi();
          data_eta = tower_geom->get_tower_geometry(keyData)->get_eta();

          sim_phi = tower_geom->get_tower_geometry(keySim)->get_phi();
          sim_eta = tower_geom->get_tower_geometry(keySim)->get_eta();
        }
        else
        {
          data_phi = tower_geomIH->get_tower_geometry(keyData)->get_phi();
          data_eta = tower_geomIH->get_tower_geometry(keyData)->get_eta();

          sim_phi = tower_geomIH->get_tower_geometry(keySim)->get_phi();
          sim_eta = tower_geomIH->get_tower_geometry(keySim)->get_eta();
        }

        if (data_phi == sim_phi && data_eta == sim_eta)
        {
          _data_towers[caloIndex]->get_tower_at_channel(channel)->set_energy(embed_E);
        }
        else
        {
          if (Verbosity())
          {
            std::cout << "eta and phi values in CEMC do not match between data and simulation, removing this event" << std::endl;
          }
          return Fun4AllReturnCodes::ABORTEVENT;
        }
      }
      else if (caloIndex == 1)
      {
        data_phi = tower_geomIH->get_tower_geometry(keyData)->get_phi();
        data_eta = tower_geomIH->get_tower_geometry(keyData)->get_eta();

        sim_phi = tower_geomIH->get_tower_geometry(keySim)->get_phi();
        sim_eta = tower_geomIH->get_tower_geometry(keySim)->get_eta();

        if (data_phi == sim_phi && data_eta == sim_eta)
        {
          _data_towers[caloIndex]->get_tower_at_channel(channel)->set_energy(embed_E);
        }
        else
        {
          if (Verbosity())
          {
            std::cout << "eta and phi values in HCALIN do not match between data and simulation, removing this event" << std::endl;
          }
          return Fun4AllReturnCodes::ABORTEVENT;
        }
      }
      else if (caloIndex == 2)
      {
        data_phi = tower_geomOH->get_tower_geometry(keyData)->get_phi();
        data_eta = tower_geomOH->get_tower_geometry(keyData)->get_eta();

        sim_phi = tower_geomOH->get_tower_geometry(keySim)->get_phi();
        sim_eta = tower_geomOH->get_tower_geometry(keySim)->get_eta();

        if (data_phi == sim_phi && data_eta == sim_eta)
        {
          _data_towers[caloIndex]->get_tower_at_channel(channel)->set_energy(embed_E);
        }
        else
        {
          if (Verbosity())
          {
            std::cout << "eta and phi values in HCALOUT do not match between data and simulation, removing this event" << std::endl;
          }
          return Fun4AllReturnCodes::ABORTEVENT;
        }
      }
    }  // end loop over channels
  }    // end loop over calorimeters

  return Fun4AllReturnCodes::EVENT_OK;
}

int caloTowerEmbed::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void caloTowerEmbed::CreateNodeTree(PHCompositeNode *topNode)
{
  std::cout << "creating node" << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *simTopNode = se->topNode("TOPSim");

  tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!tower_geom)
  {
    std::cerr << Name() << "::" << __PRETTY_FUNCTION__
              << "tower geom CEMC missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find TOWERGEOM_CEMC node");
  }

  tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!tower_geomIH)
  {
    std::cerr << Name() << "::" << __PRETTY_FUNCTION__
              << "tower geom HCALIN missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find TOWERGEOM_HCALIN node");
  }

  tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!tower_geomOH)
  {
    std::cerr << Name() << "::" << __PRETTY_FUNCTION__
              << "tower geom HCALOUT missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find TOWERGEOM_HCALOUT node");
  }

  PHNodeIterator dataIter(topNode);
  PHNodeIterator simIter(simTopNode);

  // data top node first

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(dataIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  PHCompositeNode *dstNodeSim = dynamic_cast<PHCompositeNode *>(simIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNodeSim)
  {
    std::cerr << Name() << "::" << __PRETTY_FUNCTION__
              << "DSTsim Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DSTsim node in RawTowerCalibration::CreateNodes");
  }

  for (int i = 0; i < 3; i++)
  {
    std::string detector = "CEMC";
    if (i == 0 && m_useRetower)
    {
      detector = "CEMC_RETOWER";
    }
    if (i == 1)
    {
      detector = "HCALIN";
    }
    else if (i == 2)
    {
      detector = "HCALOUT";
    }

    // data
    std::string DataTowerNodeName = "TOWERINFO_CALIB_" + detector;
    _data_towers[i] = findNode::getClass<TowerInfoContainer>(dstNode,
                                                             DataTowerNodeName);
    if (!_data_towers[i])
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << "Tower Calib Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Calib Tower node in caloTowerEmbed::CreateNodes");
    }

    // sim
    std::string SimTowerNodeName = "TOWERINFO_CALIB_" + detector;
    _sim_towers[i] = findNode::getClass<TowerInfoContainer>(dstNodeSim,
                                                            SimTowerNodeName);
    if (!_sim_towers[i])
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << "Tower Calib Sim Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Calib Tower Sim node in caloTowerEmbed::CreateNodes");
    }

  }  // end loop over calorimeter types

  return;
}

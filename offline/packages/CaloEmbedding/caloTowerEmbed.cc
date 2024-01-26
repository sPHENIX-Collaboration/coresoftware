#include "caloTowerEmbed.h"

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>

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

  if (m_dettype == CaloTowerDefs::CEMC)
    {
      m_detector = "CEMC";
    }
  else if (m_dettype == CaloTowerDefs::HCALIN)
    {
      m_detector = "HCALIN";
    }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
    {
      m_detector = "HCALOUT";
    }
  else if (m_dettype == CaloTowerDefs::ZDC)
    {
      m_detector = "ZDC";
    }
  else if (m_dettype == CaloTowerDefs::SEPD)
    {
      m_detector = "SEPD";
    }
  
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

  if (Verbosity())
    {
      std::cout << "event " << m_eventNumber << " working on " << m_detector << std::endl;
    }
  RawTowerDefs::keytype keyData = 0;
  RawTowerDefs::keytype keySim = 0;
  
  unsigned int ntowers = _data_towers->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
    {
      unsigned int data_key = _data_towers->encode_key(channel);
      unsigned int sim_key = _sim_towers->encode_key(channel);

      int ieta_data = _data_towers->getTowerEtaBin(data_key);
      int iphi_data = _data_towers->getTowerPhiBin(data_key);
      int ieta_sim = _sim_towers->getTowerEtaBin(sim_key);
      int iphi_sim = _sim_towers->getTowerPhiBin(sim_key);
      
      if (m_dettype == CaloTowerDefs::CEMC)
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
      else if (m_dettype == CaloTowerDefs::HCALIN)
	{
	  keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_data, iphi_data);
	  keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta_sim, iphi_sim);
	}
      else if (m_dettype == CaloTowerDefs::HCALOUT)
	{
	  keyData = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta_data, iphi_data);
	  keySim = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta_sim, iphi_sim);
	}
      
      TowerInfo *caloinfo_data = _data_towers->get_tower_at_channel(channel);
      TowerInfo *caloinfo_sim = _sim_towers->get_tower_at_channel(channel);

      if (m_removeBadTowers && !caloinfo_data->get_isGood())
	{
	  _data_towers->get_tower_at_channel(channel)->set_energy(0.0);
	  _data_towers->get_tower_at_channel(channel)->set_time(-11);

	  _embed_towers_out->get_tower_at_channel(channel)->set_energy(0.0);
	  _embed_towers_out->get_tower_at_channel(channel)->set_time(-11);

	  _sim_towers_out->get_tower_at_channel(channel)->set_energy(0.0);
	  _sim_towers_out->get_tower_at_channel(channel)->set_time(-11);
	  continue;
	}
      
      float data_E = caloinfo_data->get_energy();
      float sim_E = caloinfo_sim->get_energy();
      float embed_E = data_E + sim_E;
      
      float data_phi = 0.0;
      float sim_phi = 0.0;
      
      float data_eta = 0.0;
      float sim_eta = 0.0;
      
      data_phi = tower_geom->get_tower_geometry(keyData)->get_phi();
      data_eta = tower_geom->get_tower_geometry(keyData)->get_eta();
      
      sim_phi = tower_geom->get_tower_geometry(keySim)->get_phi();
      sim_eta = tower_geom->get_tower_geometry(keySim)->get_eta();
      
      if (data_phi == sim_phi && data_eta == sim_eta)
	{
	  _embed_towers_out->get_tower_at_channel(channel)->set_energy(embed_E);
	  _embed_towers_out->get_tower_at_channel(channel)->set_time(caloinfo_data->get_time());

	  _sim_towers_out->get_tower_at_channel(channel)->set_energy(sim_E);
	  _sim_towers_out->get_tower_at_channel(channel)->set_time(caloinfo_sim->get_time());
	}
      else
	{
	  if (Verbosity())
	    {
	      std::cout << "eta and phi values in " << m_detector << " do not match between data and simulation, removing this event" << std::endl;
	    }
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
           
    } // end loop over channels
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int caloTowerEmbed::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void caloTowerEmbed::CreateNodeTree(PHCompositeNode *topNode)
{
  std::cout << "creating node" << std::endl;

  std::string TowerNodeName = m_inputNodePrefix + m_detector;
  std::string GeomNodeName = "TOWERGEOM_" + m_detector;
  std::string EmbedTowerNodeName = m_inputNodePrefix + "EMBED_" + m_detector;
  std::string SimTowerNodeName = m_inputNodePrefix + "SIM_" + m_detector;
  if (m_useRetower && m_detector == "CEMC")
    {
      TowerNodeName = m_inputNodePrefix + m_detector + "_RETOWER";
      GeomNodeName = "TOWERGEOM_HCALIN";
      EmbedTowerNodeName = m_inputNodePrefix + "Embed_" + m_detector + "_RETOWER";
      SimTowerNodeName = m_inputNodePrefix + "Sim_" + m_detector + "_RETOWER";
    }

  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *simTopNode = se->topNode("TOPSim");

  tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, GeomNodeName);
  if (!tower_geom)
    {
      std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
		<< "tower geom " << GeomNodeName << " missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find " + GeomNodeName + " node");
    }

  PHNodeIterator dataIter(topNode);
  PHNodeIterator simIter(simTopNode);

  // data top node first

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(dataIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  PHCompositeNode *dstNodeSim = dynamic_cast<PHCompositeNode *>(simIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNodeSim)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DSTsim Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DSTsim node in RawTowerCalibration::CreateNodes");
  }
  
  // data
  _data_towers = findNode::getClass<TowerInfoContainer>(dstNode,TowerNodeName);
  if (!_data_towers)
    {
      std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
		<< TowerNodeName << " Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
	  "Failed to find " + TowerNodeName + " node in caloTowerEmbed::CreateNodes");
    }


  //sim
  _sim_towers = findNode::getClass<TowerInfoContainer>(dstNodeSim,TowerNodeName);
  if (!_sim_towers)
    {
      std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
		<< TowerNodeName << " Sim Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
	  "Failed to find " + TowerNodeName + " Sim node in caloTowerEmbed::CreateNodes");
      }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* caloNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", m_detector));
  
  if (!caloNode)
    {
      caloNode = new PHCompositeNode(m_detector);
      dstNode->addNode(caloNode);
    }

  TowerInfoContainer::DETECTOR towerDetector = TowerInfoContainer::DETECTOR::EMCAL;
  if (m_detector == "HCALIN" || m_detector == "HCALOUT") towerDetector = TowerInfoContainer::DETECTOR::HCAL;


  _embed_towers_out = findNode::getClass<TowerInfoContainer>(topNode, EmbedTowerNodeName);
  if (!_embed_towers_out)
    {
      _embed_towers_out = new TowerInfoContainerv2(towerDetector);
      PHIODataNode<PHObject>* embedNode =
        new PHIODataNode<PHObject>(_embed_towers_out, EmbedTowerNodeName, "PHObject");
      caloNode->addNode(embedNode);
    }


  _sim_towers_out = findNode::getClass<TowerInfoContainer>(topNode, SimTowerNodeName);
  if (!_sim_towers_out)
    {
      _sim_towers_out = new TowerInfoContainerv2(towerDetector);
      PHIODataNode<PHObject>* simNode =
        new PHIODataNode<PHObject>(_sim_towers_out, SimTowerNodeName, "PHObject");
      caloNode->addNode(simNode);
    }



  return;
}

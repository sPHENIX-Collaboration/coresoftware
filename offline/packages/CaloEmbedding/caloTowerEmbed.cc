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

  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTowerEmbed::process_event(PHCompositeNode * topNode)
{
  ++m_eventNumber;

  if (Verbosity())
  {
    std::cout << "event " << m_eventNumber << " working on " << m_detector << std::endl;
  }

  if (m_embedwaveform)
  {
    std::string ped_nodename = "PEDESTAL_" + m_detector;
    m_PedestalContainer = findNode::getClass<TowerInfoContainer>(topNode, ped_nodename);

    if (!m_PedestalContainer)
    {
      std::cout << PHWHERE << " " << ped_nodename << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }


  unsigned int ntowers = _data_towers->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {

    
    _sim_towers->get_tower_at_channel(channel)->set_status(_data_towers->get_tower_at_channel(channel)->get_status());

    TowerInfo *caloinfo_data = _data_towers->get_tower_at_channel(channel);
    TowerInfo *caloinfo_sim = _sim_towers->get_tower_at_channel(channel);

    if (m_embedwaveform)
    {
      // when the data is not ZS-ed
      if (!(caloinfo_data->get_isZS() || caloinfo_data->get_isNotInstr()))
      {
        // here we really don't want the m_samples being greater than what data has!
        for (int j = 0; j < m_nsamples; j++)
        {
          // superpose the waveforms
          caloinfo_sim->set_waveform_value(j, caloinfo_data->get_waveform_value(j) + caloinfo_sim->get_waveform_value(j));
        }
      }
      else
      {
        // if this is ZS-ed or empty(like the packet is gone)
        // we add the noise pedestal and add the post - pre to sample 6
        TowerInfo *pedestal_tower = m_PedestalContainer->get_tower_at_channel(channel);
        float pedestal_mean = 0;
        std::vector<float> m_waveform_pedestal;
        m_waveform_pedestal.resize(m_nsamples);
        // this is for pedestal scaling setting
        for (int j = 0; j < m_nsamples; j++)
        {
          m_waveform_pedestal.at(j) = (j < m_datasamples) ? pedestal_tower->get_waveform_value(j) : pedestal_tower->get_waveform_value(m_datasamples - 1);
          pedestal_mean += m_waveform_pedestal.at(j);
        }
        pedestal_mean /= m_nsamples;
        for (int j = 0; j < m_nsamples; j++)
        {
          m_waveform_pedestal.at(j) = (m_waveform_pedestal.at(j) - pedestal_mean) * m_pedestal_scale + pedestal_mean;
        }
        // add the pedestal
        for (int j = 0; j < m_nsamples; j++)
        {
          // superpose the waveforms

          caloinfo_sim->set_waveform_value(j, caloinfo_sim->get_waveform_value(j) + m_waveform_pedestal.at(j));
        }
        // add the post - pre to sample 6
        int post_pre = caloinfo_data->get_waveform_value(1) - caloinfo_data->get_waveform_value(0);

        caloinfo_sim->set_waveform_value(6, caloinfo_sim->get_waveform_value(6) + post_pre);
      }
    }
    else
    {
      float data_E = caloinfo_data->get_energy();
      float sim_E = caloinfo_sim->get_energy();
      float embed_E = data_E + sim_E;

      _sim_towers->get_tower_at_channel(channel)->set_energy(embed_E);
      _sim_towers->get_tower_at_channel(channel)->set_time(caloinfo_data->get_time());
    }

  }  // end loop over channels

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
  std::string SimTowerNodeName = TowerNodeName;
  if (m_embedwaveform)
  {
    SimTowerNodeName = m_waveformNodePrefix + m_detector;
  }
  std::string GeomNodeName = "TOWERGEOM_" + m_detector;
  if (m_useRetower && m_detector == "CEMC")
  {
    TowerNodeName = m_inputNodePrefix + m_detector + "_RETOWER";
    GeomNodeName = "TOWERGEOM_HCALIN";
  }

  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *dataTopNode = se->topNode("TOPData");

  tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, GeomNodeName);
  if (!tower_geom)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "tower geom " << GeomNodeName << " missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find " + GeomNodeName + " node");
  }

  PHNodeIterator simIter(topNode);
  PHNodeIterator dataIter(dataTopNode);

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
  _data_towers = findNode::getClass<TowerInfoContainer>(dstNode, TowerNodeName);
  if (!_data_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << TowerNodeName << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + TowerNodeName + " node in caloTowerEmbed::CreateNodes");
  }

  // sim
  _sim_towers = findNode::getClass<TowerInfoContainer>(dstNodeSim, SimTowerNodeName);
  if (!_sim_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << TowerNodeName << " Sim Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + TowerNodeName + " Sim node in caloTowerEmbed::CreateNodes");
  }

  PHNodeIterator dstIterSim(dstNodeSim);
  PHCompositeNode *caloNode = dynamic_cast<PHCompositeNode *>(dstIterSim.findFirst("PHCompositeNode", m_detector));

  if (!caloNode)
  {
    caloNode = new PHCompositeNode(m_detector);
    dstNodeSim->addNode(caloNode);
  }

  return;
}

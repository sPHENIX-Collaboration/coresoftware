#include "CaloTowerBuilder.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>

#include <caloreco/CaloWaveformProcessing.h>  // for CaloWaveformProcessing

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <iostream>  // for operator<<, endl, basic...
#include <memory>    // for allocator_traits<>::val...
#include <vector>    // for vector

//____________________________________________________________________________..
CaloTowerBuilder::CaloTowerBuilder(const std::string &name)
  : SubsysReco(name)
{
  WaveformProcessing = new CaloWaveformProcessing();
}

//____________________________________________________________________________..
CaloTowerBuilder::~CaloTowerBuilder()
{
  delete WaveformProcessing;
}

//____________________________________________________________________________..
int CaloTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  WaveformProcessing->set_softwarezerosuppression(true,40);
  if (m_dettype == CaloTowerBuilder::CEMC)
  {
    m_detector = "CEMC";
    m_packet_low = 6001;
    m_packet_high = 6128;
    // 6001, 60128
    m_nchannels = 192;
    WaveformProcessing->set_template_file("testbeam_cemc_template.root");
  }
  else if (m_dettype == CaloTowerBuilder::HCALIN)
  {
    m_packet_low = 7001;
    m_packet_high = 7008;
    m_detector = "HCALIN";
    m_nchannels = 192;
    WaveformProcessing->set_template_file("testbeam_ihcal_template.root");
  }
  else if (m_dettype == CaloTowerBuilder::HCALOUT)
  {
    m_detector = "HCALOUT";
    m_packet_low = 8001;
    m_packet_high = 8008;
    m_nchannels = 192;
    WaveformProcessing->set_template_file("testbeam_ohcal_template.root");
  }
  else if (m_dettype == CaloTowerBuilder::EPD)
  {
    m_detector = "EPD";
    m_packet_low = 9001;
    m_packet_high = 9005;
    m_nchannels = 186;
    WaveformProcessing->set_template_file("testbeam_cemc_template.root");  // place holder until we have EPD templates
  }
  WaveformProcessing->initialize_processing();
  CreateNodeTree(topNode);
  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloTowerBuilder::process_event(PHCompositeNode *topNode)
{
  std::vector<std::vector<float>> waveforms;
  if (m_isdata)
  {
    Event *_event = findNode::getClass<Event>(topNode, "PRDF");
    if (_event == nullptr)
    {
      std::cout << "CaloUnpackPRDF::Process_Event - Event not found" << std::endl;
      return -1;
    }
    if (_event->getEvtType() >= 8)  /// special event where we do not read out the calorimeters
    {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
    for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);
      if (packet)
	{
	  int nchannels = packet->iValue(0, "CHANNELS");
	  if (nchannels > m_nchannels) // packet is corrupted and reports too many channels
	    {
	      return Fun4AllReturnCodes::DISCARDEVENT;
	    }
	  for (int channel = 0; channel < nchannels; channel++)
	    {
	      std::vector<float> waveform;
	      waveform.reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{
		  waveform.push_back(packet->iValue(samp, channel));
		}
	      waveforms.push_back(waveform);
	      waveform.clear();
	    }
	  if (nchannels < m_nchannels)
	    {
	      for (int channel = 0; channel <m_nchannels - nchannels; channel++)
		{
		  std::vector<float> waveform;
		  waveform.reserve(m_nsamples);
		  for (int samp = 0; samp < m_nzerosuppsamples; samp++)
		    {
		      waveform.push_back(0);
		    }
		  waveforms.push_back(waveform);
		  waveform.clear();
		}
	    }
	  delete packet;
	}
      else // if the packet is missing treat constitutent channels as zero suppressed 
	{
	  for (int channel = 0; channel <m_nchannels; channel++) 
	    {
	        std::vector<float> waveform;
		waveform.reserve(2);
		for (int samp = 0; samp < m_nzerosuppsamples; samp++)
		  {
		    waveform.push_back(0);
		  }
		waveforms.push_back(waveform);
		waveform.clear();
	    }
	}
    }
  }
  else  // placeholder for adding simulation
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  
  std::vector<std::vector<float>> processed_waveforms = WaveformProcessing->process_waveform(waveforms);
  int n_channels = processed_waveforms.size();
  for (int i = 0; i < n_channels; i++)
  {
    m_CaloInfoContainer->get_tower_at_channel(i)->set_time(processed_waveforms.at(i).at(1));
    m_CaloInfoContainer->get_tower_at_channel(i)->set_energy(processed_waveforms.at(i).at(0));
  }
  
  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerBuilder::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);
  // DST node
  PHCompositeNode *dst_node = dynamic_cast<PHCompositeNode *>(
      nodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cout << "PHComposite node created: DST" << std::endl;
    dst_node = new PHCompositeNode("DST");
    topNode->addNode(dst_node);
  }
  // towers
  if (m_dettype == CaloTowerBuilder::CEMC)
  {
    m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::EMCAL);
  }
  else if (m_dettype == EPD)
  {
    m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::SEPD);
  }
  else
  {
    m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::HCAL);
  }

  PHIODataNode<PHObject> *emcal_towerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, "TOWERS_" + m_detector, "PHObject");
  dst_node->addNode(emcal_towerNode);
}

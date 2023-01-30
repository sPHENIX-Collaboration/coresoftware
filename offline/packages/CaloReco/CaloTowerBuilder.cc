#include "CaloTowerBuilder.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

// #include <calobase/RawTower.h>  // for RawTower
// #include <calobase/RawTowerContainer.h>
// #include <calobase/RawTowerDefs.h>  // for HCALIN, HCALOUT, CEMC

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phool/PHCompositeNode.h>
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>
// #include <TRandom3.h>

//____________________________________________________________________________..
CaloTowerBuilder::CaloTowerBuilder(const std::string &name):
 SubsysReco(name)
 , m_dettype(CaloTowerBuilder::CEMC)
 , m_CaloInfoContainer(0)
 , m_detector("CEMC")
 , m_packet_low(6017)
 , m_packet_high(6032)
 , m_nsamples(16)
 , m_isdata(true)
{
  std::cout << "CaloTowerBuilder::CaloTowerBuilder(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloTowerBuilder::~CaloTowerBuilder()
{
  std::cout << "CaloTowerBuilder::~CaloTowerBuilder() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  // rnd = new TRandom3();

  WaveformProcessing = new CaloWaveformProcessing();
  WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);

  if (m_dettype == CaloTowerBuilder::CEMC)
    {
      m_detector = "CEMC";
      m_packet_low = 6017;
      m_packet_high = 6032;
      // 6001, 60128
      WaveformProcessing->set_template_file("testbeam_cemc_template.root");
    }
  else if (m_dettype == CaloTowerBuilder::HCALIN)
    {
      m_packet_low  = 7001;
      m_packet_high = 7032;
      m_detector = "HCALIN";
      WaveformProcessing->set_template_file("testbeam_ihcal_template.root");
   }
  else if (m_dettype == CaloTowerBuilder::HCALOUT)
    {
      m_detector = "HCALOUT";
      m_packet_low = 8001;
      m_packet_high = 8032;
      WaveformProcessing->set_template_file("testbeam_ohcal_template.root");
   }
  else if (m_dettype == CaloTowerBuilder::EPD)
    {
      m_detector = "EPD";
      m_packet_low = 9001;
      m_packet_high = 9016;
      WaveformProcessing->set_template_file("testbeam_cemc_template.root"); // place holder until we have EPD templates
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
      if (_event == 0)
      	{
      	  std::cout << "CaloUnpackPRDF::Process_Event - Event not found" << std::endl;
      	  return -1;
      	}
      if ( _event->getEvtType() >= 8)/// special event where we do not read out the calorimeters
      	{
      	  return Fun4AllReturnCodes::DISCARDEVENT;
      	}
      for ( int pid = m_packet_low; pid <= m_packet_high; pid++)
      	{ 
      	  Packet *packet = _event->getPacket(pid);
      	  if (!packet)
      	    {
      	      return Fun4AllReturnCodes::DISCARDEVENT;
      	    }      
      	  for ( int channel = 0; channel <  packet->iValue(0,"CHANNELS"); channel++)
      	    {
      	      std::vector<float> waveform;
      	      for (int samp = 0; samp < m_nsamples;samp++)
      	      	{
      	      	  waveform.push_back(packet->iValue(samp,channel));	      
      	      	}
      	      waveforms.push_back(waveform);
      	      waveform.clear();
      	    }
	  delete packet;
      	}
    }
  else // placeholder for adding simulation 
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

  std::vector<std::vector<float>> processed_waveforms =  WaveformProcessing->process_waveform(waveforms);
  int n_channels = processed_waveforms.size();
  for (int i = 0 ; i < n_channels;i++)
    {
      m_CaloInfoContainer->at(i)->set_time(processed_waveforms.at(i).at(1));
      m_CaloInfoContainer->at(i)->set_energy(processed_waveforms.at(i).at(0));
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
      m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainerv1::DETECTOR::EMCAL);
    }
  else if (m_dettype == EPD)
    {
      m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainerv1::DETECTOR::SEPD);
    }
  else
    {
      m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainerv1::DETECTOR::HCAL);  
    }
  
  PHIODataNode<PHObject> *emcal_towerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, Form("TOWERS_%s",m_detector.c_str()), "PHObject");
  // PHIODataNode<PHObject> *emcal_towerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, "TOWERS_"+m_detector, "PHObject");
  dst_node->addNode(emcal_towerNode);
}







#include "LL1PacketGetter.h"

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
LL1PacketGetter::LL1PacketGetter(const std::string &name, const std::string &trigger, const std::string &ll1)
  : SubsysReco(name)
  , m_trigger(trigger)
  , m_ll1(ll1)
  , m_ll1out(nullptr)
  , m_packet_low(INT_MIN)
  , m_packet_high(INT_MIN)
  , m_nsamples(20)
  , m_nchannels(60)
  , m_isdata(true)
{

  m_prim_map[TriggerDefs::DetectorId::noneDId] = 0;
  m_prim_map[TriggerDefs::DetectorId::emcalDId] = 384;
  m_prim_map[TriggerDefs::DetectorId::hcalinDId] = 24;
  m_prim_map[TriggerDefs::DetectorId::hcaloutDId] = 24;
  m_prim_map[TriggerDefs::DetectorId::mbdDId] = 4;

  _verbose = 0;
}

//____________________________________________________________________________..
LL1PacketGetter::~LL1PacketGetter()
{

}


//____________________________________________________________________________..
int LL1PacketGetter::InitRun(PHCompositeNode *topNode)
{

  m_triggerid = TriggerDefs::GetTriggerId(m_trigger);
  m_detectorid = TriggerDefs::GetDetectorId(m_ll1);
  m_primitiveid = TriggerDefs::GetPrimitiveId(m_ll1);

  m_nprimitives = m_prim_map[m_detectorid];
  
  switch (m_detectorid) {
  case TriggerDefs::DetectorId::hcalDId:
  case TriggerDefs::DetectorId::hcalinDId:
  case TriggerDefs::DetectorId::hcaloutDId:
      std::cout <<"HCAL Packet Getter."<<std::endl;
      m_packet_low = 14003;
      m_packet_high = 14003;
      m_nchannels_per_primitive = 16;
      m_ntriggerwords = 1;
      break;
  case TriggerDefs::DetectorId::mbdDId:
      std::cout <<"MBD Packet Getter"<<std::endl;
      m_packet_low = 14002;
      m_packet_high = 14002;
      m_nchannels_per_primitive = 13;
      m_ntriggerwords = 8;
      break;

  default:
      std::cout <<"This Packet Getter not implemented yet."<<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;

      break;

    }
  m_nchannels = m_nchannels_per_primitive*m_nprimitives + m_ntriggerwords;
  
  if (CreateNodeTree(topNode)) return Fun4AllReturnCodes::ABORTRUN;

  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LL1PacketGetter::process_event(PHCompositeNode *topNode)
{
  //  std::cout << m_detector_type <<std::endl;
  switch(m_detectorid){
  case TriggerDefs::DetectorId::hcalinDId :
  case TriggerDefs::DetectorId::hcaloutDId :
  case TriggerDefs::DetectorId::hcalDId :
    {
      
      m_ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_RAW_HCAL");
      if (!m_ll1out)
	{
	  std::cout << "HCAL LL1 data not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }
  case TriggerDefs::DetectorId::mbdDId :
    {
      
      m_ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_RAW_MBD");
      if (!m_ll1out)
	{
	  std::cout << "MBD LL1 data not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }
  default:
      std::cout <<"This Packet Getter not implemented yet."<<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;

      break;

  }
  TriggerDefs::TriggerSumKey sumkey = 0;
  TriggerDefs::TriggerPrimKey primkey;
  _trigger_primitives = m_ll1out->GetTriggerPrimitiveContainer();
  if (m_isdata)
  {
    Event *_event = findNode::getClass<Event>(topNode, "PRDF");
    if (_event == nullptr)
    {
      std::cout << "LL1UnpackPRDF::Process_Event - Event not found" << std::endl;
      return -1;
    }
    if (_event->getEvtType() >= 8)  /// special event where we do not read out the calorimeters
    {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
    unsigned int clk;
    unsigned int evt;
    
    for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);
      if (packet)
	{
	  evt = packet->iValue(0, "EVTNR");
	  clk = packet->iValue(0, "CLOCK");
	  m_ll1out->set_event_number(evt);
	  m_ll1out->set_clock_number(clk);

	  int nchannels = packet->iValue(0, "CHANNELS");
	  if (nchannels > m_nchannels) // packet is corrupted and reports too many channels
	    {
	      return Fun4AllReturnCodes::DISCARDEVENT;
	    }

	  for (int iprim = 0; iprim <m_nprimitives ; iprim++)
	    {

	      if (m_detectorid == TriggerDefs::DetectorId::hcalDId)
		{
		  primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId((iprim < 12 ? "HCALOUT": "HCALIN")), TriggerDefs::GetPrimitiveId(m_ll1), iprim);
		}
	      else 
		{
		  primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), TriggerDefs::GetPrimitiveId(m_ll1), iprim);
		}	
	      _trigger_primitive = new TriggerPrimitive(primkey);
	      for (int channel = 0; channel < m_nchannels_per_primitive; channel++)
		    {
		      sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), TriggerDefs::GetPrimitiveId(m_ll1), iprim, channel);
		      _sum = new std::vector<unsigned int>();
		      _sum->reserve(m_nsamples);
		      for (int samp = 0; samp < m_nsamples; samp++)
			{
			  _sum->push_back(static_cast<unsigned int>(packet->iValue(samp, iprim*m_nchannels_per_primitive+channel)));
			}

		      _trigger_primitive->add_sum(sumkey, _sum);
		    }
	      _trigger_primitives->add_primitive(primkey, _trigger_primitive);
	    }
	  for (int channel = 0; channel < m_ntriggerwords;channel++)
	    {
	      _sum = new std::vector<unsigned int>();
	      _sum->reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{
		  _sum->push_back(static_cast<unsigned int>(packet->iValue(samp, m_nprimitives*m_nchannels_per_primitive+channel)));
		}

	      m_ll1out->add_word(channel, _sum);
	      
	    }
	
	  delete packet;
	}
      else // if the packet is missing treat constitutent channels as zero suppressed 
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }
  }
  else  // placeholder for adding simulation
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (_verbose) m_ll1out->identify();

  return Fun4AllReturnCodes::EVENT_OK;
}

int LL1PacketGetter::CreateNodeTree(PHCompositeNode *topNode)
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
  PHNodeIterator dstIter(dst_node);

  PHCompositeNode *ll1Node = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "LL1"));
  if (!ll1Node)
    {
      ll1Node = new PHCompositeNode("LL1");
      dst_node->addNode(ll1Node);
    }
  switch (m_detectorid) 
    {
    case TriggerDefs::DetectorId::mbdDId :
      {
	LL1Outv2 *ll1out = findNode::getClass<LL1Outv2>(ll1Node, "LL1OUT_RAW_MBD");
	if (!ll1out)
	  {
	    ll1out = new LL1Outv2(m_trigger, m_ll1);
	    PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, "LL1OUT_RAW_MBD", "PHObject");
	    ll1Node->addNode(LL1OutNode);
	  }
	break;
      }
    case TriggerDefs::DetectorId::hcalDId :
    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
      {
	LL1Outv2 *ll1out = findNode::getClass<LL1Outv2>(ll1Node, "LL1OUT_RAW_HCAL");
	if (!ll1out)
	  {
	    ll1out = new LL1Outv2(m_trigger, m_ll1);
	    PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, "LL1OUT_RAW_HCAL", "PHObject");
	    ll1Node->addNode(LL1OutNode);
	  }
	break;
      }
    default:
      std::cout <<"This Packet Getter not implemented yet."<<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
      
      break;

    }
  return Fun4AllReturnCodes::EVENT_OK;
}

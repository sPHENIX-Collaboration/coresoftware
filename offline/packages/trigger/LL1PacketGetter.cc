#include "LL1PacketGetter.h"
#include "LL1Outv1.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

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
#include <limits>
#include <memory>  // for allocator_traits<>::val...
#include <vector>  // for vector

//____________________________________________________________________________..
LL1PacketGetter::LL1PacketGetter(const std::string &name, const std::string &trigger, const std::string &ll1)
  : SubsysReco(name)
  , m_trigger(trigger)
  , m_ll1(ll1)
  , m_ll1out(nullptr)
  , m_packet_low(std::numeric_limits<int>::min())
  , m_packet_high(std::numeric_limits<int>::min())
  , m_nchannels(60)
  , m_isdata(true)
  , m_no_ll1out(false)
{
  m_triggerid = TriggerDefs::TriggerId::noneTId;
  m_primitiveid = TriggerDefs::PrimitiveId::nonePId;
  m_detectorid = TriggerDefs::DetectorId::noneDId;

  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(0, 0);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalinDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcaloutDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcalDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcalinDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcaloutDId)] = std::make_pair(13002, 13002);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(13010, 13025);
  m_packet_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::mbdTId, TriggerDefs::DetectorId::mbdDId)] = std::make_pair(13002, 13002);

  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(0, 0);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalinDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcaloutDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::noneDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcalDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcalinDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::hcaloutDId)] = std::make_pair(16, 24);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::emcalDId)] = std::make_pair(24, 16);
  m_prim_sum_map[TriggerDefs::getTriggerKey(TriggerDefs::TriggerId::mbdTId, TriggerDefs::DetectorId::mbdDId)] = std::make_pair(4, 13);

  m_word_map[TriggerDefs::TriggerId::noneTId] = 0;
  m_word_map[TriggerDefs::TriggerId::jetTId] = 12 * 32;
  m_word_map[TriggerDefs::TriggerId::mbdTId] = 8;
  m_word_map[TriggerDefs::TriggerId::pairTId] = 0;
  m_word_map[TriggerDefs::TriggerId::photonTId] = 12 * 32;
}

//____________________________________________________________________________..
int LL1PacketGetter::InitRun(PHCompositeNode *topNode)
{
  m_triggerid = TriggerDefs::GetTriggerId(m_trigger);
  m_detectorid = TriggerDefs::GetDetectorId(m_ll1);
  m_triggerkey = TriggerDefs::getTriggerKey(m_triggerid, m_detectorid);
  m_primitiveid = TriggerDefs::GetPrimitiveId(m_ll1);

  if (strcmp(m_trigger.c_str(), "NONE") == 0)
  {
    m_triggerprimitive_nodename = "TRIGGERPRIMITIVES_RAW_" + m_ll1;
    m_triggerprimitive_ll1_nodename = "TRIGGERPRIMITIVES_RAW_" + m_ll1 + "_LL1";
    m_no_ll1out = true;
  }
  else
  {
    m_ll1_nodename = "LL1OUT_RAW_" + m_trigger;
    m_triggerprimitive_nodename = "TRIGGERPRIMITIVES_RAW_" + m_trigger;
  }

  m_packet_low = m_packet_map[m_triggerkey].first;
  m_packet_high = m_packet_map[m_triggerkey].second;
  m_nprimitives = m_prim_sum_map[m_triggerkey].first;
  m_nchannels_per_primitive = m_prim_sum_map[m_triggerkey].second;
  m_ntriggerwords = m_word_map[m_triggerid];

  m_nchannels = m_nchannels_per_primitive * m_nprimitives + m_ntriggerwords;

  if (CreateNodeTree(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LL1PacketGetter::process_event(PHCompositeNode *topNode)
{
  m_trigger_primitives = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_triggerprimitive_nodename);
  if (!m_trigger_primitives)
  {
    std::cout << m_triggerprimitive_nodename << " not found - Fatal Error" << std::endl;
    exit(1);
  }
  if (!m_no_ll1out)
  {
    m_trigger_primitives->setTriggerType(m_triggerid);

    m_ll1out = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);
    if (!m_ll1out)
    {
      std::cout << m_ll1_nodename << " not found - Fatal Error" << std::endl;
      exit(1);
    }
  }
  else
  {
    m_trigger_primitives->setTriggerType(m_triggerid);
    m_trigger_primitives_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_triggerprimitive_ll1_nodename);
    if (!m_trigger_primitives_ll1)
    {
      std::cout << m_triggerprimitive_ll1_nodename << " not found - Fatal Error" << std::endl;
      exit(1);
    }
  }
  TriggerDefs::TriggerSumKey sumkey = 0;
  TriggerDefs::TriggerPrimKey primkey;

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
    int nsamples;
    //int monitor;
    for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);

      if (packet)
      {
        evt = packet->iValue(0, "EVTNR");
        clk = packet->iValue(0, "CLOCK");
        nsamples = packet->iValue(0, "SAMPLES");
        //monitor = packet->iValue(0, "MONITOR");

        if (!m_no_ll1out)
        {
          m_ll1out->set_event_number(evt);
          m_ll1out->set_clock_number(clk);
        }

        int nchannels = packet->iValue(0, "CHANNELS");
        int ntriggerwords = packet->iValue(0, "TRIGGERWORDS");

        if (nchannels > m_nchannels)  // packet is corrupted and reports too many channels
        {
          return Fun4AllReturnCodes::DISCARDEVENT;
        }

	TriggerDefs::PrimitiveId primid = TriggerDefs::GetPrimitiveId(m_trigger);
	if (primid == TriggerDefs::PrimitiveId::nonePId)
	{
	  primid = TriggerDefs::GetPrimitiveId(m_ll1);
	}

        for (int iprim = 0; iprim < m_nprimitives; iprim++)
        {

	  
          primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), primid, m_nprimitives * (pid - m_packet_low) + iprim);

          _trigger_primitive = new TriggerPrimitivev1(primkey);

          for (int channel = 0; channel < m_nchannels_per_primitive; channel++)
          {
	    sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), primid, m_nprimitives * (pid - m_packet_low) + iprim, channel);
	    
            _sum = new std::vector<unsigned int>();
            _sum->reserve(nsamples);
            for (int samp = 0; samp < nsamples; samp++)
            {
              if (pid == 13002)
	      {
		_sum->push_back(static_cast<unsigned int>(packet->iValue(samp, iprim*2 + channel/12 + 32*(channel%12))));
	      }
	      else
	      {
		_sum->push_back(static_cast<unsigned int>(packet->iValue(samp, iprim * m_nchannels_per_primitive + channel)));
	      }
            }

            _trigger_primitive->add_sum(sumkey, _sum);
          }
          m_trigger_primitives->add_primitive(primkey, _trigger_primitive);
        }

        if (!m_no_ll1out)
        {
          for (int channel = 0; channel < ntriggerwords; channel++)
          {
            _sum = new std::vector<unsigned int>();
            _sum->reserve(nsamples);
            for (int samp = 0; samp < nsamples; samp++)
            {
              _sum->push_back(static_cast<unsigned int>(packet->iValue(samp, m_nprimitives * m_nchannels_per_primitive + channel)));
            }

            m_ll1out->add_word(((unsigned int) (channel % 32) & 0xffffU) + (((unsigned int) (channel / 32) & 0xffffU) << 16U), _sum);
          }
        }
        else
        {
          TriggerDefs::TriggerId tid = TriggerDefs::TriggerId::jetTId;
          TriggerDefs::PrimitiveId prid = TriggerDefs::PrimitiveId::jetPId;
          TriggerDefs::DetectorId did = TriggerDefs::DetectorId::emcalDId;

	  
          m_trigger_primitives_ll1->setTriggerType(tid);
          primkey = TriggerDefs::getTriggerPrimKey(tid, did, prid, pid - m_packet_low);

          _trigger_primitive = new TriggerPrimitivev1(primkey);
	  
          for (int channel = 0; channel < 24; channel++)
          {
            sumkey = TriggerDefs::getTriggerSumKey(tid, did, prid, pid - m_packet_low, channel);
            _sum = new std::vector<unsigned int>();
            _sum->reserve(nsamples);

            for (int samp = 0; samp < nsamples; samp++)
            {
              _sum->push_back(static_cast<unsigned int>(packet->iValue(samp, m_nprimitives * m_nchannels_per_primitive + channel)));
            }
            _trigger_primitive->add_sum(sumkey, _sum);
          }
          m_trigger_primitives_ll1->add_primitive(primkey, _trigger_primitive);
        }

        delete packet;
      }
      // else // if the packet is missing treat constitutent channels as zero suppressed
      //   {
      //     for (int iprim = 0; iprim < m_nprimitives; iprim++)
      // 	{
      // 	  if (m_detectorid == TriggerDefs::DetectorId::hcalDId)
      // 	    {
      // 	      primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId((iprim < 12 ? "HCALOUT": "HCALIN")), TriggerDefs::GetPrimitiveId(m_ll1), m_nprimitives*(pid - m_packet_low) + iprim);
      // 	    }
      // 	  else
      // 	    {
      // 	      primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), TriggerDefs::GetPrimitiveId(m_ll1), m_nprimitives*(pid - m_packet_low) + iprim);
      // 	    }

      // 	  _trigger_primitive = new TriggerPrimitivev1(primkey);

      // 	  for (int channel = 0; channel < m_nchannels_per_primitive; channel++)
      // 	    {
      // 	      sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId(m_ll1), TriggerDefs::GetPrimitiveId(m_ll1), m_nprimitives*(pid - m_packet_low) + iprim, channel);
      // 	      _sum = new std::vector<unsigned int>();
      // 	      _sum->reserve(nsamples);
      // 	      for (int samp = 0; samp < nsamples; samp++)
      // 		{
      // 		  _sum->push_back(0); //static_cast<unsigned int>(packet->iValue(samp, iprim*m_nchannels_per_primitive+channel)));
      // 		}

      // 	      _trigger_primitive->add_sum(sumkey, _sum);
      // 	    }
      // 	  m_trigger_primitives->add_primitive(primkey, _trigger_primitive);
      // 	}

      //     for (int channel = 0; channel < m_ntriggerwords;channel++)
      // 	{
      // 	  _sum = new std::vector<unsigned int>();
      // 	  _sum->reserve(nsamples);
      // 	  for (int samp = 0; samp < nsamples; samp++)
      // 	    {
      // 	      _sum->push_back(0);
      // 	    }

      // 	  m_ll1out->add_word(unsigned int)(channel), _sum);

      // 	}

      //   }
    }
  }

  if (Verbosity())
  {
    if (m_no_ll1out)
    {
      std::cout << " PRIMITIVES LL1: " << std::endl;
      m_trigger_primitives_ll1->identify();
    }
    else
    {
      std::cout << " LL1 OUT: " << std::endl;
      m_ll1out->identify();
    }
  }
  if (Verbosity())
  {
    std::cout << "PRIMITIVES: " << std::endl;
    m_trigger_primitives->identify();
  }

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

  if (!m_no_ll1out)
  {
    LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
    if (!ll1out)
    {
      ll1out = new LL1Outv1(m_trigger, m_ll1);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  else
  {
    TriggerPrimitiveContainer *primoutll1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_triggerprimitive_ll1_nodename);
    if (!primoutll1)
    {
      primoutll1 = new TriggerPrimitiveContainerv1();
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(primoutll1, m_triggerprimitive_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  TriggerPrimitiveContainer *primout = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_triggerprimitive_nodename);
  if (!primout)
  {
    primout = new TriggerPrimitiveContainerv1();
    PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(primout, m_triggerprimitive_nodename, "PHObject");
    ll1Node->addNode(LL1OutNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

#include "SingleLL1TriggerInput.h"

#include "Fun4AllPrdfInputTriggerManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/LL1PacketContainerv1.h>
#include <ffarawobjects/LL1Packetv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/packet.h>  // for Packet

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream<...
#include <iterator>  // for reverse_iterator
#include <limits>    // for numeric_limits
#include <memory>
#include <set>
#include <utility>  // for pair

// it is 8 packets for the ll1, this number needs to be npackets+1
// so it doesn't trigger the warning and exit. Setting it to 10
static const int NLL1PACKETS = 19;

SingleLL1TriggerInput::SingleLL1TriggerInput(const std::string &name)
  : SingleTriggerInput(name)
{
  SubsystemEnum(InputManagerType::LL1);
  plist = new Packet *[NLL1PACKETS];  // eight packets for the ll1 in each file
}

SingleLL1TriggerInput::~SingleLL1TriggerInput()
{
  CleanupUsedPackets(std::numeric_limits<int>::max());
  // some events are already in the m_EventStack but they haven't been put
  // into the m_PacketMap
  while (m_EventStack.begin() != m_EventStack.end())
  {
    m_EventStack.erase(m_EventStack.begin());
  }
  delete[] plist;
}

void SingleLL1TriggerInput::FillPool(const unsigned int keep)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }
  while (GetSomeMoreEvents(keep))
  {
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt.reset(GetEventiterator()->getNextEvent());
    }
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, NLL1PACKETS);
    if (npackets >= NLL1PACKETS)
    {
      std::cout << PHWHERE << " Packets array size " << NLL1PACKETS
                << " too small for " << Name()
                << ", increase NLL1PACKETS and rebuild" << std::endl;
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      int packet_id = plist[i]->getIdentifier();
      // The call to  EventNumberOffset(identifier) will initialize it to our default (zero) if it wasn't set already
      // if we encounter a misalignemt, the Fun4AllPrdfInputTriggerManager will adjust this. But the event
      // number of the adjustment depends on its pooldepth. Events in its pools will be moved to the correct slots
      // and only when the pool gets refilled, this correction kicks in
      // SO DO NOT BE CONFUSED when printing this out - seeing different events where this kicks in
      int CorrectedEventSequence = EventSequence + EventNumberOffset(packet_id);
      if (Verbosity() > 2)
      {
        plist[i]->identify();
      }

      // by default use previous bco clock for gtm bco
      LL1Packet *newhit = new LL1Packetv1();
      int nr_channels = plist[i]->iValue(0, "CHANNELS");
      int nr_samples = plist[i]->iValue(0, "SAMPLES");
      uint64_t gtm_bco = plist[i]->iValue(0, "CLOCK");
      // offline packet content
      newhit->setEvtSequence(CorrectedEventSequence);
      newhit->setIdentifier(packet_id);
      newhit->setHitFormat(plist[i]->getHitFormat());
      newhit->setBCO(gtm_bco);
      newhit->setPacketEvtSequence(plist[i]->iValue(0, "EVTNR"));
      // ll1 packet additions
      newhit->setNrSamples(nr_samples);
      newhit->setNrChannels(nr_channels);
      newhit->setTriggerWords(plist[i]->iValue(0, "TRIGGERWORDS"));
      newhit->setSlotNr(plist[i]->iValue(0, "SLOTNR"));
      newhit->setCardNr(plist[i]->iValue(0, "CARDNR"));
      newhit->setMonitor(plist[i]->iValue(0, "MONITOR"));
      newhit->setFemWords(plist[i]->iValue(0, "FEMWORDS"));
      newhit->setFibers(plist[i]->iValue(0, "FIBERS"));
      newhit->setSums(plist[i]->iValue(0, "SUMS"));
      for (int ichan = 0; ichan < nr_channels; ichan++)
      {
        for (int isamp = 0; isamp < nr_samples; isamp++)
        {
          if (isamp >= newhit->getMaxNumSamples() || ichan >= newhit->getMaxNumChannels())
          {
            std::cout << "Packet: " << newhit->getIdentifier()
                      << ", samples: " << isamp
                      << ", channels: " << ichan << std::endl;
            gSystem->Exit(1);
          }
          else
          {
            newhit->setSample(ichan, isamp, plist[i]->iValue(isamp, ichan));
          }
        }
      }
      // newhit->identify();
      //       newhit->dump();
      if (Verbosity() > 2)
      {
        std::cout << PHWHERE << "corrected evtno: " << CorrectedEventSequence
                  << ", original evtno: " << EventSequence
                  << ", bco: 0x" << std::hex << gtm_bco << std::dec
                  << std::endl;
      }
      if (TriggerInputManager())
      {
        TriggerInputManager()->AddLL1Packet(CorrectedEventSequence, newhit);
      }
      m_PacketMap[CorrectedEventSequence].push_back(newhit);
      m_EventStack.insert(CorrectedEventSequence);
      if (ddump_enabled())
      {
        ddumppacket(plist[i]);
      }
      delete plist[i];
    }
  }
}

void SingleLL1TriggerInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_PacketMap)
    {
      std::cout << PHWHERE << "Event: " << bcliter.first << std::endl;
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_EventStack)
    {
      std::cout << PHWHERE << "stacked event: " << iter << std::endl;
    }
  }
}

void SingleLL1TriggerInput::CleanupUsedPackets(const int eventno)
{
  std::vector<int> toclearevents;
  for (const auto &iter : m_PacketMap)
  {
    if (iter.first <= eventno)
    {
      for (auto pktiter : iter.second)
      {
        delete pktiter;
      }
      toclearevents.push_back(iter.first);
    }
    else
    {
      break;
    }
  }

  for (auto iter : toclearevents)
  {
    m_EventStack.erase(iter);
    m_PacketMap.erase(iter);
  }
}

void SingleLL1TriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  int currentevent = *m_EventStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentevent);
  return;
}

void SingleLL1TriggerInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "LL1"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("LL1");
    dstNode->addNode(detNode);
  }
  LL1PacketContainer *ll1packetcont = findNode::getClass<LL1PacketContainer>(detNode, "LL1Packets");
  if (!ll1packetcont)
  {
    ll1packetcont = new LL1PacketContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(ll1packetcont, "LL1Packets", "PHObject");
    detNode->addNode(newNode);
  }
}

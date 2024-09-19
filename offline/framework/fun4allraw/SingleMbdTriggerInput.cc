#include "SingleMbdTriggerInput.h"

#include "Fun4AllPrdfInputTriggerManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

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

SingleMbdTriggerInput::SingleMbdTriggerInput(const std::string &name)
  : SingleTriggerInput(name)
{
  SubsystemEnum(InputManagerType::MBD);
}

SingleMbdTriggerInput::~SingleMbdTriggerInput()
{
  CleanupUsedPackets(std::numeric_limits<int>::max());
  // some events are already in the m_EventStack but they haven't been put
  // into the m_PacketMap
  while (m_EventStack.begin() != m_EventStack.end())
  {
    m_EventStack.erase(m_EventStack.begin());
  }
}

void SingleMbdTriggerInput::FillPool(const unsigned int keep)
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
    if (EventSequence < SkipToEvent())
    {
      continue;
    }
    std::vector<Packet *> pktvec = evt->getPacketVector();
    for (auto packet : pktvec)
    {
      int packet_id = packet->getIdentifier();
      // The call to  EventNumberOffset(identifier) will initialize it to our default if it wasn't set already
      // if we encounter a misalignemt, the Fun4AllPrdfInputTriggerManager will adjust this. But the event
      // number of the adjustment depends on its pooldepth. Events in its pools will be moved to the correct slots
      // and only when the pool gets refilled, this correction kicks in
      // SO DO NOT BE CONFUSED when printing this out - seeing different events where this kicks in
      int CorrectedEventSequence = EventSequence + EventNumberOffset(packet_id);
      if (Verbosity() > 2)
      {
        packet->identify();
      }

      CaloPacket *newhit = new CaloPacketv1();
      int nr_modules = packet->iValue(0, "NRMODULES");
      int nr_channels = packet->iValue(0, "CHANNELS");
      int nr_samples = packet->iValue(0, "SAMPLES");
      if (nr_modules > 3)
      {
        std::cout << PHWHERE << " too many modules, need to adjust arrays" << std::endl;
        gSystem->Exit(1);
      }

      uint64_t gtm_bco = packet->lValue(0, "CLOCK");
      newhit->setNrModules(nr_modules);
      newhit->setNrSamples(nr_samples);
      newhit->setNrChannels(nr_channels);
      newhit->setBCO(gtm_bco);
      newhit->setPacketEvtSequence(packet->iValue(0, "EVTNR"));
      newhit->setIdentifier(packet_id);
      newhit->setHitFormat(packet->getHitFormat());
      newhit->setEvtSequence(CorrectedEventSequence);
      newhit->setEvenChecksum(packet->iValue(0, "EVENCHECKSUM"));
      newhit->setCalcEvenChecksum(packet->iValue(0, "CALCEVENCHECKSUM"));
      newhit->setOddChecksum(packet->iValue(0, "ODDCHECKSUM"));
      newhit->setCalcOddChecksum(packet->iValue(0, "CALCODDCHECKSUM"));
      newhit->setModuleAddress(packet->iValue(0, "MODULEADDRESS"));
      newhit->setDetId(packet->iValue(0, "DETID"));
      for (int ifem = 0; ifem < nr_modules; ifem++)
      {
        newhit->setFemClock(ifem, packet->iValue(ifem, "FEMCLOCK"));
        newhit->setFemEvtSequence(ifem, packet->iValue(ifem, "FEMEVTNR"));
        newhit->setFemSlot(ifem, packet->iValue(ifem, "FEMSLOT"));
        newhit->setChecksumLsb(ifem, packet->iValue(ifem, "CHECKSUMLSB"));
        newhit->setChecksumMsb(ifem, packet->iValue(ifem, "CHECKSUMMSB"));
        newhit->setCalcChecksumLsb(ifem, packet->iValue(ifem, "CALCCHECKSUMLSB"));
        newhit->setCalcChecksumMsb(ifem, packet->iValue(ifem, "CALCCHECKSUMMSB"));
      }
      for (int ipmt = 0; ipmt < nr_channels; ipmt++)
      {
        // store pre/post only for suppressed channels, the array in the packet routines is not
        // initialized so reading pre/post for not zero suppressed channels returns garbage
        bool isSuppressed = packet->iValue(ipmt, "SUPPRESSED");
        newhit->setSuppressed(ipmt, isSuppressed);
        if (isSuppressed)
        {
          newhit->setPre(ipmt, packet->iValue(ipmt, "PRE"));
          newhit->setPost(ipmt, packet->iValue(ipmt, "POST"));
        }
        else
        {
          for (int isamp = 0; isamp < nr_samples; isamp++)
          {
            newhit->setSample(ipmt, isamp, packet->iValue(isamp, ipmt));
          }
        }
      }
      if (Verbosity() > 2)
      {
        std::cout << PHWHERE << "corrected evtno: " << CorrectedEventSequence
                  << ", original evtno: " << EventSequence
                  << ", bco: 0x" << std::hex << gtm_bco << std::dec
                  << std::endl;
      }
      if (TriggerInputManager())
      {
        TriggerInputManager()->AddMbdPacket(CorrectedEventSequence, newhit);
      }
      m_PacketMap[CorrectedEventSequence].push_back(newhit);
      m_EventStack.insert(CorrectedEventSequence);
      if (ddump_enabled())
      {
        ddumppacket(packet);
      }
      delete packet;
    }
  }
}

void SingleMbdTriggerInput::Print(const std::string &what) const
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

void SingleMbdTriggerInput::CleanupUsedPackets(const int eventno)
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

void SingleMbdTriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  int currentevent = *m_EventStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentevent);
  return;
}

void SingleMbdTriggerInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "MBD"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("MBD");
    dstNode->addNode(detNode);
  }
  CaloPacketContainer *mbdpacketcont = findNode::getClass<CaloPacketContainer>(detNode, "MBDPackets");
  if (!mbdpacketcont)
  {
    mbdpacketcont = new CaloPacketContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(mbdpacketcont, "MBDPackets", "PHObject");
    detNode->addNode(newNode);
  }
}

#include "SingleHcalTriggerInput.h"

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

// it is 8 packets for the hcal, this number needs to be npackets+1
// so it doesn't trigger the warning and exit. Setting it to 10
static const int NHCALPACKETS = 10;

SingleHcalTriggerInput::SingleHcalTriggerInput(const std::string &name)
  : SingleTriggerInput(name)
{
  SubsystemEnum(InputManagerType::HCAL);
  plist = new Packet *[NHCALPACKETS];  // eight packets for the hcal in each file
  LocalPoolDepth(3);
}

SingleHcalTriggerInput::~SingleHcalTriggerInput()
{
  CleanupUsedLocalPackets(std::numeric_limits<int>::max());
  CleanupUsedPackets(std::numeric_limits<int>::max());
  // some events are already in the m_EventStack but they haven't been put
  // into the m_PacketMap
  while (m_EventStack.begin() != m_EventStack.end())
  {
    m_EventStack.erase(m_EventStack.begin());
  }
  delete[] plist;
}

void SingleHcalTriggerInput::FillPool(const unsigned int keep)
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
    int npackets = evt->getPacketList(plist, NHCALPACKETS);
    if (npackets >= NHCALPACKETS)
    {
      std::cout << PHWHERE << " Packets array size " << NHCALPACKETS
                << " too small for " << Name()
                << ", increase NHCALPACKETS and rebuild" << std::endl;
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
      CaloPacket *newhit = new CaloPacketv1();
      int nr_modules = plist[i]->iValue(0, "NRMODULES");
      int nr_channels = plist[i]->iValue(0, "CHANNELS");
      int nr_samples = plist[i]->iValue(0, "SAMPLES");
      if (nr_modules > newhit->getMaxNumModules())
      {
        std::cout << PHWHERE << " too many modules " << nr_modules << ", max is "
                  << newhit->getMaxNumModules() << ", need to adjust arrays" << std::endl;
        gSystem->Exit(1);
      }
      if (nr_channels > newhit->getMaxNumChannels())
      {
        std::cout << PHWHERE << " too many channels " << nr_channels << ", max is "
                  << newhit->getMaxNumChannels() << ", need to adjust arrays" << std::endl;
        gSystem->Exit(1);
      }
      if (nr_samples > newhit->getMaxNumSamples())
      {
        std::cout << PHWHERE << " too many samples " << nr_samples << ", max is "
                  << newhit->getMaxNumSamples() << ", need to adjust arrays" << std::endl;
        gSystem->Exit(1);
      }
      uint64_t gtm_bco = plist[i]->lValue(0, "CLOCK");
      newhit->setNrModules(nr_modules);
      newhit->setNrSamples(nr_samples);
      newhit->setNrChannels(nr_channels);
      newhit->setBCO(gtm_bco);
      newhit->setPacketEvtSequence(plist[i]->iValue(0, "EVTNR"));
      newhit->setIdentifier(packet_id);
      newhit->setHitFormat(plist[i]->getHitFormat());
      newhit->setEvtSequence(EventSequence);
      newhit->setEvenChecksum(plist[i]->iValue(0, "EVENCHECKSUM"));
      newhit->setCalcEvenChecksum(plist[i]->iValue(0, "CALCEVENCHECKSUM"));
      newhit->setOddChecksum(plist[i]->iValue(0, "ODDCHECKSUM"));
      newhit->setCalcOddChecksum(plist[i]->iValue(0, "CALCODDCHECKSUM"));
      newhit->setModuleAddress(plist[i]->iValue(0, "MODULEADDRESS"));
      newhit->setDetId(plist[i]->iValue(0, "DETID"));
      for (int ifem = 0; ifem < nr_modules; ifem++)
      {
        newhit->setFemClock(ifem, plist[i]->iValue(ifem, "FEMCLOCK"));
        newhit->setFemEvtSequence(ifem, plist[i]->iValue(ifem, "FEMEVTNR"));
        newhit->setFemSlot(ifem, plist[i]->iValue(ifem, "FEMSLOT"));
        newhit->setChecksumLsb(ifem, plist[i]->iValue(ifem, "CHECKSUMLSB"));
        newhit->setChecksumMsb(ifem, plist[i]->iValue(ifem, "CHECKSUMMSB"));
        newhit->setCalcChecksumLsb(ifem, plist[i]->iValue(ifem, "CALCCHECKSUMLSB"));
        newhit->setCalcChecksumMsb(ifem, plist[i]->iValue(ifem, "CALCCHECKSUMMSB"));
      }
      for (int ipmt = 0; ipmt < nr_channels; ipmt++)
      {
        // store pre/post only for suppressed channels, the array in the packet routines is not
        // initialized so reading pre/post for not zero suppressed channels returns garbage
        bool isSuppressed = plist[i]->iValue(ipmt, "SUPPRESSED");
        newhit->setSuppressed(ipmt, isSuppressed);
        if (isSuppressed)
        {
          newhit->setPre(ipmt, plist[i]->iValue(ipmt, "PRE"));
          newhit->setPost(ipmt, plist[i]->iValue(ipmt, "POST"));
        }
        else
        {
          for (int isamp = 0; isamp < nr_samples; isamp++)
          {
            newhit->setSample(ipmt, isamp, plist[i]->iValue(isamp, ipmt));
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
      m_LocalPacketMap[CorrectedEventSequence].push_back(newhit);
      m_EventStack.insert(CorrectedEventSequence);
      if (ddump_enabled())
      {
        ddumppacket(plist[i]);
      }
      delete plist[i];
    }
    if (m_LocalPacketMap.size() >= LocalPoolDepth())
    {
      CheckFEMEventNumber();
    }
    while (m_LocalPacketMap.size() > LocalPoolDepth())
    {
      std::set<int> events;
      auto nh = m_LocalPacketMap.begin()->second;
      //      std::cout << "Pushing event " << m_LocalPacketMap.begin()->first << " from local packet map to packet map" << std::endl;
      events.insert(m_LocalPacketMap.begin()->first);
      m_PacketMap[m_LocalPacketMap.begin()->first] = std::move(nh);
      m_LocalPacketMap.erase(m_LocalPacketMap.begin());
      // copy reference bco if we have a FEM problem
      if (FEMClockProblemFlag())
      {
        std::vector<OfflinePacket *> badpackets;
        uint64_t refbco = std::numeric_limits<uint64_t>::max();
        uint64_t fallbackrefbco = std::numeric_limits<uint64_t>::max();
        for (const auto &iter : m_PacketMap)
        {
          if (events.find(iter.first) == events.end())
          {
            //	    std::cout << "event " << iter.first << " already bco treated" << std::endl;
            continue;
          }
          for (auto pktiter : iter.second)
          {
            if (pktiter->getIdentifier() == ClockReferencePacket())
            {
              refbco = pktiter->getBCO();
            }
            else if (m_BadBCOPacketSet.find(pktiter->getIdentifier()) == m_BadBCOPacketSet.end())
            {
              fallbackrefbco = pktiter->getBCO();  // all bcos are identical so we can pcik any for fallback
            }
            else
            {
              //	      std::cout << "found bad packet " << pktiter->getIdentifier() << std::endl;
              badpackets.push_back(pktiter);
            }
          }
          if (refbco == std::numeric_limits<uint64_t>::max())
          {
            static int count = 0;
            if (count < 1000)
            {
              std::cout << PHWHERE << ": crap that didn't work, could not locate reference packet" << std::endl;
              count++;
            }
            refbco = fallbackrefbco;
          }
          for (auto pktiter : badpackets)
          {
            if (Verbosity() > 2)
            {
              std::cout << "event " << iter.first << " Setting packet " << pktiter->getIdentifier() << " to bco " << std::hex
                        << refbco << std::dec << std::endl;
            }
            pktiter->setBCO(refbco);
          }
        }
      }
    }
    //    Print("PACKETMAP");
    if (TriggerInputManager())
    {
      for (const auto &evtiter : m_PacketMap)
      {
        for (auto pktiter : evtiter.second)
        {
          CaloPacket *calpacket = dynamic_cast<CaloPacket *>(pktiter);
          if (calpacket)
          {
            //	    std::cout << "pushing event " << evtiter.first << " to combiner" << std::endl;
            TriggerInputManager()->AddHcalPacket(evtiter.first, calpacket);
          }
          else
          {
            static int count = 0;
            if (count < 1000)
            {
              std::cout << PHWHERE << " dynamic cast from offline to calo packet failed??? here is its identify():" << std::endl;
              count++;
            }
            pktiter->identify();
          }
        }
      }
    }
  }
}

void SingleHcalTriggerInput::Print(const std::string &what) const
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
  if (what == "LOCALMAP")
  {
    for (auto &iter : m_LocalPacketMap)
    {
      std::cout << "LocalMap Event " << iter.first << std::endl;
      for (auto pktiter : iter.second)
      {
        std::cout << "Packet " << pktiter->getIdentifier()
                  << ", BCO: " << std::hex << pktiter->getBCO() << std::dec
                  << ", FEM: " << std::hex << pktiter->iValue(0, "FEMCLOCK") << std::dec << std::endl;
      }
    }
  }
  if (what == "PACKETMAP")
  {
    for (auto &iter : m_PacketMap)
    {
      std::cout << "PacketMap Event " << iter.first << std::endl;
      for (auto pktiter : iter.second)
      {
        std::cout << "Packet " << pktiter->getIdentifier()
                  << ", BCO: " << std::hex << pktiter->getBCO() << std::dec
                  << ", FEM: " << std::hex << pktiter->iValue(0, "FEMCLOCK") << std::dec << std::endl;
      }
    }
  }
}

void SingleHcalTriggerInput::CleanupUsedLocalPackets(const int eventno)
{
  std::vector<int> toclearevents;
  for (const auto &iter : m_LocalPacketMap)
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
    // std::set::erase returns number of elements deleted if the key does not exist, it just returns 0
    m_EventStack.erase(iter);
    m_LocalPacketMap.erase(iter);
  }
}

void SingleHcalTriggerInput::CleanupUsedPackets(const int eventno)
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
    // std::set::erase returns number of elements deleted if the key does not exist, it just returns 0
    m_EventStack.erase(iter);
    m_PacketMap.erase(iter);
  }
}

void SingleHcalTriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  int currentevent = *m_EventStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentevent);
  return;
}

void SingleHcalTriggerInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "HCAL"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("HCAL");
    dstNode->addNode(detNode);
  }
  CaloPacketContainer *hcalpacketcont = findNode::getClass<CaloPacketContainer>(detNode, "HCALPackets");
  if (!hcalpacketcont)
  {
    hcalpacketcont = new CaloPacketContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hcalpacketcont, "HCALPackets", "PHObject");
    detNode->addNode(newNode);
  }
}
void SingleHcalTriggerInput::CheckFEMClock()
{
  // lets check in the first event if this is actually needed
  auto first_event = m_LocalPacketMap.begin();
  //    std::cout << "Event " << first_event->first << std::endl;
  std::map<uint64_t, std::set<int>> pktbcomap;
  uint64_t ref_femclk = std::numeric_limits<uint64_t>::max();
  std::map<uint64_t, unsigned int> bcocount;
  for (auto pktiter : first_event->second)
  {
    //      std::cout << "packet id: " << pktiter->getIdentifier() << " size: " <<  first_event->second.size() << std::endl;
    std::set<uint64_t> femclockset;
    for (int i = 0; i < pktiter->iValue(0, "NRMODULES"); i++)
    {
      uint64_t femclk = pktiter->iValue(i, "FEMCLOCK");
      bcocount[femclk]++;
      if (Verbosity() > 21)
      {
        std::cout << "packet id: " << pktiter->getIdentifier() << " packet clock: 0x" << std::hex << pktiter->iValue(0, "CLOCK")
                  << " FEMClock: 0x" << femclk << std::dec << std::endl;
      }
      femclockset.insert(femclk);
      if (ref_femclk == std::numeric_limits<uint64_t>::max())
      {
        ref_femclk = femclk;
      }
      else
      {
        if (ref_femclk != femclk)
        {
          if (Verbosity() > 1)
          {
            std::cout << "Event " << first_event->first << " FEM Clock mismatch for packet " << pktiter->getIdentifier() << std::endl;
            std::cout << "ref fem clk: 0x" << std::hex << ref_femclk << ", femclk: 0x"
                      << femclk << std::dec << std::endl;
          }
        }
      }
      if (femclockset.size() > 1)
      {
        static int count = 0;
        if (count < 1000)
        {
          std::cout << PHWHERE << " FEM Clocks differ for packet " << pktiter->getIdentifier()
                    << ", found " << femclockset.size() << " different ones" << std::endl;
          for (auto &iter : femclockset)
          {
            std::cout << "0x" << std::hex << iter << std::dec << std::endl;
          }
          count++;
        }
      }
    }
    pktbcomap[*femclockset.begin()].insert(pktiter->getIdentifier());
  }
  //    std::cout << "Map size " << bcocount.size() << std::endl;
  if (bcocount.size() < 2)
  {
    //      std::cout << "all good" << std::endl;
    return;
  }
  static int count = 0;
  if (count < 1000)
  {
    std::cout << PHWHERE << " FEM clocks are off, found " << bcocount.size() << " different ones, here we go ..." << std::endl;
    count++;
  }
  // set our FEM Clock Problem flag, since we need to copy clocks in the Fill loop
  SetFEMClockProblemFlag();
  // find good bco (which will give us the haystack) and bad packets
  if (Verbosity() > 1)
  {
    std::cout << "LocalPacketMap size: " << m_LocalPacketMap.size()
              << ", pool depth: " << LocalPoolDepth() << std::endl;
  }
  if (m_LocalPacketMap.size() < LocalPoolDepth())
  {
    // std::cout << "cache is not deep enough, this should never happen, size of local packet map: "
    // 		<< m_LocalPacketMap.size() << ", pool depth: " << LocalPoolDepth() << std::endl;
    return;
  }
  // first find the reference bco (majority of packets until I get the Master from JeaBeom
  uint64_t goodfembco = std::numeric_limits<uint64_t>::max();
  unsigned int maxnumpackets = 0;
  for (auto bcoiter : bcocount)
  {
    if (bcoiter.second > maxnumpackets)
    {
      maxnumpackets = bcoiter.second;
      goodfembco = bcoiter.first;
    }
    // 	std::cout << "bco 0x" << std::hex << bcoiter.first << " shows up " << std::dec << bcoiter.second << std::endl;
  }
  int refpacketid = *pktbcomap.find(goodfembco)->second.begin();
  if (Verbosity() > 1)
  {
    std::cout << "Use packet " << refpacketid << " for reference bco 0x"
              << std::hex << goodfembco << std::dec << std::endl;
  }
  SetClockReferencePacket(refpacketid);
  pktbcomap.erase(goodfembco);
  for (const auto &badpktmapiter : pktbcomap)
  {
    for (auto badpktiter : badpktmapiter.second)
    {
      //	   std::cout << "bad packet " << badpktiter << std::endl;
      if (TriggerInputManager())
      {
        TriggerInputManager()->AddFEMProblemPacket(badpktiter);
      }
      m_BadBCOPacketSet.insert(badpktiter);
    }
  }
  std::vector<uint64_t> HayStack;
  std::map<int, std::vector<uint64_t>> NeedleMap;
  m_EventRefBCO.clear();  // this is used if we need to reshuffle the BCO
  for (auto &iter : m_LocalPacketMap)
  {
    //    std::cout << "handling Event " << iter->first << std::endl;
    for (auto pktiter : iter.second)
    {
      if (pktiter->getIdentifier() == refpacketid)
      {
        // just pick the first one, we have already checked that they are identical
        uint64_t femclk = pktiter->iValue(0, "FEMCLOCK");
        HayStack.push_back(femclk);
        m_EventRefBCO[iter.first] = pktiter->getBCO();
      }
      else if (m_BadBCOPacketSet.find(pktiter->getIdentifier()) != m_BadBCOPacketSet.end())
      {
        uint64_t femclk = pktiter->iValue(0, "FEMCLOCK");
        NeedleMap[pktiter->getIdentifier()].push_back(femclk);
      }
    }
  }
  if (Verbosity() > 1)
  {
    for (auto bco : HayStack)
    {
      std::cout << "Haystack : 0x" << std::hex << bco << std::dec << std::endl;
    }
  }
  for (const auto &needleiter : NeedleMap)
  {
    std::vector needle = needleiter.second;
    needle.pop_back();
    if (Verbosity() > 1)
    {
      std::cout << "Packet " << needleiter.first << std::endl;
      for (auto bco : needle)
      {
        std::cout << "Needle: 0x" << std::hex << bco << std::dec << std::endl;
      }
    }
    auto it = std::search(HayStack.begin(), HayStack.end(), needle.begin(), needle.end());
    if (it != HayStack.end())  // found the needle in the haystack at offset position
    {
      int position = std::distance(HayStack.begin(), it);
      //     std::cout << "found needle at " << position << std::endl;
      AdjustEventNumberOffset(needleiter.first, position);
      ShiftEvents(needleiter.first, position);
    }
  }
  return;
}

int SingleHcalTriggerInput::ShiftEvents(int pktid, int offset)
{
  std::vector<int> eventnumbers;
  for (auto pktmapiter = m_LocalPacketMap.rbegin(); pktmapiter != m_LocalPacketMap.rend(); ++pktmapiter)
  {
    eventnumbers.push_back(pktmapiter->first);
  }
  for (auto evtnumiter : eventnumbers)
  {
    auto &pktmapiter = m_LocalPacketMap[evtnumiter];

    int newevent = evtnumiter + offset;
    for (unsigned int i = 0; i < pktmapiter.size(); ++i)
    {
      auto packet = pktmapiter[i];
      if (packet->getIdentifier() == pktid)
      {
        if (Verbosity() > 1)
        {
          std::cout << "moving packet " << packet->getIdentifier() << " from position " << i
                    << " from event " << evtnumiter << " to event " << newevent << std::endl;
        }
        auto bcotmpiter = m_EventRefBCO.find(newevent);
        if (bcotmpiter != m_EventRefBCO.end())
        {
          packet->setBCO(bcotmpiter->second);
        }
        else
        {
          packet->setBCO(std::numeric_limits<uint64_t>::max());
        }

        m_LocalPacketMap[newevent].push_back(packet);
        pktmapiter.erase(pktmapiter.begin() + i);
        break;
      }
    }
    if (Verbosity() > 1)
    {
      for (auto iter : m_LocalPacketMap[evtnumiter])
      {
        std::cout << "local packetmap after erase: " << iter->getIdentifier() << std::endl;
      }
    }
  }
  //  Print("LOCALMAP");
  return 0;
}

void SingleHcalTriggerInput::CheckFEMEventNumber()
{
  // lets check in the first event if this is actually needed
  auto first_event = m_LocalPacketMap.begin();
  //    std::cout << "Event " << first_event->first << std::endl;
  std::map<int, std::set<int>> pktevtnummap;
  int ref_femevtnum = std::numeric_limits<int>::max();
  std::map<int, unsigned int> evtnumcount;
  for (auto pktiter : first_event->second)
  {
    //      std::cout << "packet id: " << pktiter->getIdentifier() << " size: " <<  first_event->second.size() << std::endl;
    std::set<int> femevtnumset;
    for (int i = 0; i < pktiter->iValue(0, "NRMODULES"); i++)
    {
      int femevtnum = pktiter->iValue(i, "FEMEVTNR");
      evtnumcount[femevtnum]++;
      if (Verbosity() > 21)
      {
        std::cout << "packet id: " << pktiter->getIdentifier() << " packet clock: 0x" << std::hex << pktiter->iValue(0, "CLOCK")
                  << " FEM EvtNo: " << std::dec << femevtnum << std::endl;
      }
      femevtnumset.insert(femevtnum);
      if (ref_femevtnum == std::numeric_limits<int>::max())
      {
        ref_femevtnum = femevtnum;
      }
      else
      {
        if (ref_femevtnum != femevtnum)
        {
          if (Verbosity() > 1)
          {
            std::cout << "Event " << first_event->first << " FEM Event Number  mismatch for packet " << pktiter->getIdentifier() << std::endl;
            std::cout << "ref fem evt: " << ref_femevtnum << ", femevtnum: "
                      << femevtnum << std::endl;
          }
        }
      }
      if (femevtnumset.size() > 1)
      {
        static int count = 0;
        if (count < 1000)
        {
          std::cout << PHWHERE << " FEM Event Numbers differ for packet " << pktiter->getIdentifier()
                    << ", found " << femevtnumset.size() << " different ones" << std::endl;
          for (auto &iter : femevtnumset)
          {
            std::cout << iter << std::endl;
          }
          count++;
        }
      }
    }
    pktevtnummap[*femevtnumset.begin()].insert(pktiter->getIdentifier());
  }
  //    std::cout << "Map size " << evtnumcount.size() << std::endl;
  if (evtnumcount.size() < 2)
  {
    //      std::cout << "all good" << std::endl;
    return;
  }
  static int count = 0;
  if (count < 1000)
  {
    std::cout << PHWHERE << " FEM clocks are off, found " << evtnumcount.size() << " different ones, here we go ..." << std::endl;
    count++;
  }
  // set our FEM Clock Problem flag, since we need to copy clocks in the Fill loop
  SetFEMClockProblemFlag();
  // find good bco (which will give us the haystack) and bad packets
  if (Verbosity() > 1)
  {
    std::cout << "LocalPacketMap size: " << m_LocalPacketMap.size()
              << ", pool depth: " << LocalPoolDepth() << std::endl;
  }
  if (m_LocalPacketMap.size() < LocalPoolDepth())
  {
    // std::cout << "cache is not deep enough, this should never happen, size of local packet map: "
    // 		<< m_LocalPacketMap.size() << ", pool depth: " << LocalPoolDepth() << std::endl;
    return;
  }
  // first find the reference bco (majority of packets until I get the Master from JeaBeom
  int goodfemevtnum = std::numeric_limits<int>::max();
  unsigned int maxnumpackets = 0;
  for (auto bcoiter : evtnumcount)
  {
    if (bcoiter.second > maxnumpackets)
    {
      maxnumpackets = bcoiter.second;
      goodfemevtnum = bcoiter.first;
    }
    // 	std::cout << "bco 0x" << std::hex << bcoiter.first << " shows up " << std::dec << bcoiter.second << std::endl;
  }
  int refpacketid = *pktevtnummap.find(goodfemevtnum)->second.begin();
  if (Verbosity() > 1)
  {
    std::cout << "Use packet " << refpacketid << " for reference bco 0x"
              << std::hex << goodfemevtnum << std::dec << std::endl;
  }
  SetClockReferencePacket(refpacketid);
  pktevtnummap.erase(goodfemevtnum);
  for (const auto &badpktmapiter : pktevtnummap)
  {
    for (auto badpktiter : badpktmapiter.second)
    {
      //	   std::cout << "bad packet " << badpktiter << std::endl;
      if (TriggerInputManager())
      {
        TriggerInputManager()->AddFEMProblemPacket(badpktiter);
      }
      m_BadBCOPacketSet.insert(badpktiter);
    }
  }
  std::vector<int> HayStack;
  std::map<int, std::vector<int>> NeedleMap;
  m_EventRefBCO.clear();  // this is used if we need to reshuffle the BCO
  for (auto &iter : m_LocalPacketMap)
  {
    //    std::cout << "handling Event " << iter->first << std::endl;
    for (auto pktiter : iter.second)
    {
      if (pktiter->getIdentifier() == refpacketid)
      {
        // just pick the first one, we have already checked that they are identical
        int femevtnum = pktiter->iValue(0, "FEMEVTNR");
        HayStack.push_back(femevtnum);
        m_EventRefBCO[iter.first] = pktiter->getBCO();
      }
      else if (m_BadBCOPacketSet.find(pktiter->getIdentifier()) != m_BadBCOPacketSet.end())
      {
        int femevtnum = pktiter->iValue(0, "FEMEVTNR");
        NeedleMap[pktiter->getIdentifier()].push_back(femevtnum);
      }
    }
  }
  if (Verbosity() > 1)
  {
    for (auto evtno : HayStack)
    {
      std::cout << "Haystack : " << evtno << std::endl;
    }
  }
  for (const auto &needleiter : NeedleMap)
  {
    std::vector needle = needleiter.second;
    needle.pop_back();
    if (Verbosity() > 1)
    {
      std::cout << "Packet " << needleiter.first << std::endl;
      for (auto evtno : needle)
      {
        std::cout << "Needle: " << evtno << std::endl;
      }
    }
    auto it = std::search(HayStack.begin(), HayStack.end(), needle.begin(), needle.end());
    if (it != HayStack.end())  // found the needle in the haystack at offset position
    {
      int position = std::distance(HayStack.begin(), it);
      //     std::cout << "found needle at " << position << std::endl;
      AdjustEventNumberOffset(needleiter.first, position);
      ShiftEvents(needleiter.first, position);
    }
  }
  return;
}

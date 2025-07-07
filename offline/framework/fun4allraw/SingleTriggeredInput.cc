#include "SingleTriggeredInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <ranges>
#include <set>
#include <utility>  // for pair
#include <vector>

SingleTriggeredInput::SingleTriggeredInput(const std::string &name)
  : Fun4AllBase(name)
{
  m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
}

SingleTriggeredInput::~SingleTriggeredInput()
{
  while (m_EventDeque.begin() != m_EventDeque.end())
  {
    delete m_EventDeque.front();
    m_EventDeque.pop_front();
  }
  delete m_EventIterator;
}

int SingleTriggeredInput::fileopen(const std::string &filenam)
{
  std::cout << PHWHERE << "trying to open " << filenam << std::endl;
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << FileName()
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  std::string fname = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << FileName() << std::endl;
  }
  int status = 0;
  m_EventIterator = new fileEventiterator(fname.c_str(), status);
  m_EventsThisFile = 0;
  if (status)
  {
    delete m_EventIterator;
    m_EventIterator = nullptr;
    std::cout << PHWHERE << Name() << ": could not open file " << fname << std::endl;
    return -1;
  }
  IsOpen(1);
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int SingleTriggeredInput::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  delete m_EventIterator;
  m_EventIterator = nullptr;
  IsOpen(0);
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  UpdateFileList();
  return 0;
}

int SingleTriggeredInput::FillEventVector(int index)
{
  while (GetEventIterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return -1;
    }
  }
  //  std::cout << "index: " << index << std::endl;
  //  if (!m_EventDeque.empty())
  // if (!NeedsRefill())
  // {
  //   return 0;
  // }
  int i = index;
  // uint64_t tmp = m_bclkarray[pooldepth];
  // m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  // m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
  // m_bclkarray[0] = tmp;
  // std::cout << std::hex << "m_bclkarray[0]: 0x" << m_bclkarray[0] << ", ind " << std::dec << pooldepth
  //  	    << std::hex << ", 0x" << m_bclkarray[pooldepth] << std::dec << std::endl;
  size_t loopcount{0};
  while (loopcount < 1)
  {
    Event *evt{nullptr};
    do
    {
      evt = GetEventIterator()->getNextEvent();
      while (!evt)
      {
        fileclose();
        if (!OpenNextFile())
        {
          FilesDone(1);
          if (Verbosity() > 0)
          {
            std::cout << "no more events to read, deque depth: " << m_EventDeque.size() << std::endl;
          }
          return -1;
        }
        evt = GetEventIterator()->getNextEvent();
      }
      m_Gl1_SkipEvents--;
      // in run2pp up to run 48706 (first ok physics run with our seletcion) the gl1 packet numbers
      // are off (each by 2?) which triggers the gl1 skip event
      // This patch calculates the current clock diffand compares to the gl1. If they match - the gl1
      // must not be skipped. Chances of a random match are low - this would require a 32 bit rollover which
      // at 10MHz is > 400 sec between events and you have to hit the exact bunch crossing
      // even if this happens - the next even will just fail the clock check
      if (m_Gl1_SkipEvents > 0 && Gl1Input() != this)
      {
        uint64_t tmpclkdiff = calccurrentclockdiff(index, evt);
        if (tmpclkdiff == Gl1Input()->get_clkdiff(index))
        {
          if (Verbosity() > 0)
          {
            std::cout << Name() << " clock diff works out, not skipping " << evt->getEvtSequence() << std::endl;
          }
          m_Gl1_SkipEvents = 0;
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << Name() << " Skipping Event " << evt->getEvtSequence() << std::endl;
          }
        }
      }
    } while (m_Gl1_SkipEvents > 0);
    //    std::cout << "donefillig: " << DoneFilling() << std::endl;
    m_EventsThisFile++;
    // std::cout << Name() << ": ";
    // evt->identify();
    if (evt->getEvtType() != DATAEVENT)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << " dropping non data event: " << evt->getEvtSequence() << std::endl;
      }
      delete evt;
      continue;
    }
    // this is for tests
    if (m_ProblemEvent >= 0)
    {
      if (m_ProblemEvent == 0)
      {
        m_ProblemEvent--;
        delete evt;
        continue;
      }

      m_ProblemEvent--;
    }
    //    std::cout << Name() << std::endl;
    //    evt->identify();
    evt->convert();
    if (firstcall)
    {
      CreateDSTNodes(evt);
      firstcall = false;
    }
    //    std::cout << Name() << ":filling index " << i+1 << " with " << std::hex << bla << std::dec << std::endl;
    uint64_t myClock = GetClock(evt);
    if (Verbosity() > 0)
    {
      std::cout << Name() << " evt# " << evt->getEvtSequence() << " clock: 0x" << std::hex
                << myClock << std::dec << std::endl;
    }
    if (myClock == std::numeric_limits<uint64_t>::max())
    {
      std::cout << Name() << " dropping bad clock event: " << evt->getEvtSequence() << std::endl;
      delete evt;
      continue;
    }
    m_bclkarray[i + 1] = myClock;
    if (m_Event == 0)
    {
      m_bclkdiffarray[i] = 0;
    }
    else
    {
      if (m_bclkarray[i + 1] < m_bclkarray[i])
      {
        m_bclkdiffarray[i] = m_bclkarray[i + 1] + 0x100000000 - m_bclkarray[i];
      }
      else
      {
        m_bclkdiffarray[i] = m_bclkarray[i + 1] - m_bclkarray[i];
      }
      // this is just a safeguard addressing the previous wrong handling of the rollover
      // I leave it here just in case this happens again
      if (m_bclkdiffarray[i] != (m_bclkdiffarray[i] & 0xFFFFFFFF))
      {
        std::cout << Name() << " This should not happen: Found upper 32bits set: 0x" << std::hex << m_bclkdiffarray[i]
                  << " for event # " << std::dec << evt->getEvtSequence() << std::endl;
        std::cout << std::hex << "current clk: 0x" << myClock << " corrected: 0x" << m_bclkarray[i + 1]
                  << ", previous clock : 0x" << m_bclkarray[i] << std::dec << std::endl;
        m_bclkdiffarray[i] &= 0xFFFFFFFF;
      }
    }
    // std::cout << Name() << " evt# " << evt->getEvtSequence() << std::hex
    // 	      << " prev clk: 0x" << m_bclkarray[i]
    // 	      << " clock: 0x" << m_bclkarray[i + 1]
    // 	      << " clock diff: 0x" << m_bclkdiffarray[i]
    // 	      << std::dec << std::endl;

    m_Event++;
    // std::cout << "m_bclkarray[" << i<< "]: 0x" << std::hex << m_bclkarray[i] << std::dec
    // 	      << ", m_bclkarray[" << i << "+1]: 0x" << std::hex << m_bclkarray[i+1] << std::dec
    // 	      << ", m_bclkdiffarray[" << i << "]: 0x" << std::hex << m_bclkdiffarray[i] << std::dec << std::endl;
    m_EventDeque.push_back(evt);

    loopcount++;
  }
  // std::cout << Name() << std::endl;
  // for (auto iter : m_bclkdiffarray)
  // {
  //   std::cout << "0x" << std::hex << iter << std::dec << std::endl;
  // }
  return m_EventDeque.size();
}

uint64_t SingleTriggeredInput::GetClock(Event *evt)
{
  int refclockset = false;
  uint64_t clock = std::numeric_limits<uint64_t>::max();
  std::vector<Packet *> pktvec = evt->getPacketVector();

  size_t offset = 1;
  for (auto *iter : pktvec | std::views::drop(offset))
  {
    if (Verbosity() > 1)
    {
      std::cout << "checking packet " << iter->getIdentifier() << std::endl;
    }
    uint64_t pktclock = static_cast<uint64_t>(iter->lValue(0, "CLOCK") & 0xFFFFFFFF);  // NOLINT (hicpp-signed-bitwise)
    if (iter->getStatus())
    {
      std::cout << "Event " << evt->getEvtSequence() << " Packet " << iter->getIdentifier()
                << " has non zero status: " << iter->getStatus()
                << " dropping it from clock check" << std::endl;
      continue;
    }
    if (!refclockset)
    {
      clock = pktclock;
      refclockset = true;
      continue;
    }
    if (clock != pktclock)  // NOLINT (hicpp-signed-bitwise)
    {
      static int icnt = 0;
      if (icnt < 100)
      {
        std::cout << "clock problem for packet " << iter->getIdentifier() << std::hex
                  << " clock value: 0x" << pktclock << std::dec
                  << " ref clock: 0x" << clock << std::dec
                  << ", status: " << iter->getStatus() << std::endl;
        icnt++;
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "Getting clock from packet " << pktvec[0]->getIdentifier() << std::endl;
  }
  //     uint64_t clock = pktvec[0]->lValue(0, "CLOCK");
  //     pktvec[0]->identify();
  // std::cout << "pkt event: " << pktvec[0]->iValue(0, "EVTNR") << ", clock: 0x"
  // 	       << std::hex << pktvec[0]->lValue(0, "CLOCK") << std::dec << std::endl;
  for (auto *iter : pktvec)
  {
    delete iter;
  }
  return clock;
}

void SingleTriggeredInput::FillPool(int index)
{
  if (AllDone() || EventAlignmentProblem())  // no more files and all events read or alignment problem
  {
    return;
  }
  if (!FilesDone())
  {
    //    std::cout << "need to skip " << Gl1Input()->SkipEvents() << std::endl;
    m_Gl1_SkipEvents = Gl1Input()->SkipEvents();
    int iret = FillEventVector(index);
    if (Verbosity() > 0)
    {
      std::cout << Name() << " return FillEventVector( " << index << "): " << iret << std::endl;
    }
    //     // if (iret != 0)
    //     // {
    //       // for (const auto *itertst = clkdiffbegin(); itertst != clkdiffend(); ++itertst)
    //       // {
    //       // 	std::cout << std::hex << "blkdiff: 0x" << itertst << std::dec << std::endl;
    //       // }
    // //RunCheck();
    //    }
  }
  return;
}

void SingleTriggeredInput::CreateDSTNodes(Event *evt)
{
  std::string CompositeNodeName = "Packets";
  if (KeepMyPackets())
  {
    CompositeNodeName = "PacketsKeep";
  }
  PHNodeIterator iter(m_topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    m_topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", CompositeNodeName));
  if (!detNode)
  {
    detNode = new PHCompositeNode(CompositeNodeName);
    dstNode->addNode(detNode);
  }
  std::vector<Packet *> pktvec = evt->getPacketVector();
  for (auto *piter : pktvec)
  {
    int packet_id = piter->getIdentifier();
    m_PacketSet.insert(packet_id);
    std::string PacketNodeName = std::to_string(packet_id);
    CaloPacket *calopacket = findNode::getClass<CaloPacket>(detNode, PacketNodeName);
    if (!calopacket)
    {
      calopacket = new CaloPacketv1();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(calopacket, PacketNodeName, "PHObject");
      detNode->addNode(newNode);
    }
    delete piter;
  }
}

int SingleTriggeredInput::FemEventNrClockCheck(OfflinePacket *pkt)
{
  CaloPacket *calopkt = dynamic_cast<CaloPacket *>(pkt);
  if (!calopkt)  // this is a dumb check, should really segfault
  {
    std::cout << PHWHERE << ": packet null pointer, returning zero" << std::endl;
    return 0;
  }
  if (calopkt->getStatus())
  {
    return -1;
  }
  // make sure all clocks of the FEM are fine,
  int nrModules = calopkt->iValue(0, "NRMODULES");
  std::set<int> EventNoSet;
  for (int j = 0; j < nrModules; j++)
  {
    EventNoSet.insert(calopkt->iValue(j, "FEMEVTNR"));
  }
  size_t femeventnumbers = EventNoSet.size();
  if (femeventnumbers > 1)
  {
    int goodfemevent = 0;  // store the good event number so we can insert it in the set, but only if the clock counters agree
    if (femeventnumbers == 2)
    {
      // find the outlier if we have a 2:1 decision
      std::map<int, int> EventMap;
      std::map<int, int> BadModuleMap;
      for (int j = 0; j < nrModules; j++)
      {
        EventMap[calopkt->iValue(j, "FEMEVTNR")]++;
        BadModuleMap[calopkt->iValue(j, "FEMEVTNR")] = j;
      }
      for (const auto iter : EventMap)
      {
        if (iter.second == 1)
        {
          calopkt->setFemStatus(BadModuleMap[iter.first], CaloPacket::BAD_EVENTNR);
        }
        else
        {
          goodfemevent = iter.first;
        }
      }
    }
    else  // all event numbers are different - mark all bad
    {
      for (int j = 0; j < nrModules; j++)
      {
        calopkt->setFemStatus(j, CaloPacket::BAD_EVENTNR);
      }
    }
    // at least one packet (6024) has a stuck bit in the fem event nr, check fem clock counter in this case
    // if they are identical FEM is good (not checked if the FEM clock is stuck though)
    std::set<int> FemClockSet;
    for (int j = 0; j < nrModules; j++)
    {
      FemClockSet.insert(calopkt->iValue(j, "FEMCLOCK"));
    }
    if (FemClockSet.size() == 1)
    {
      static int icnt = 0;
      if (icnt < 10)
      {
        icnt++;
        std::cout << "Packet " << calopkt->getIdentifier() << " has not unique event numbers"
                  << " but FEM Clock counters are identical" << std::endl;
      }
      if (goodfemevent > 0)
      {
        m_FEMEventNrSet.insert(goodfemevent);
      }
      return 1;
    }
    // now lets find which one is the outlier
    static int icnt = 0;
    if (icnt < 1000)
    {
      icnt++;
      std::cout << "resetting packet " << calopkt->getIdentifier()
                << " with fem event and clock mismatch" << std::endl;
      std::map<int, int> ClockMap;
      std::map<int, int> EventMap;
      for (int j = 0; j < nrModules; j++)
      {
        EventMap[calopkt->iValue(j, "FEMEVTNR")]++;
        ClockMap[calopkt->iValue(j, "FEMCLOCK")]++;
      }
      for (const auto iterA : EventMap)
      {
        std::cout << "Event Nr : " << iterA.first << " shows up " << iterA.second << " times"
                  << std::hex << ", Event Nr 0x" << iterA.first << std::dec << std::endl;
      }
      for (const auto iterA : ClockMap)
      {
        std::cout << "Clock : 0x" << std::hex << iterA.first << std::dec
                  << " shows up " << iterA.second << " times" << std::endl;
      }
    }
    return -1;
  }
  m_FEMEventNrSet.insert(*(EventNoSet.begin()));
  return 0;
}

void SingleTriggeredInput::dumpdeque()
{
  const auto *iter1 = clkdiffbegin();
  const auto *iter2 = Gl1Input()->clkdiffbegin();
  while (iter1 != clkdiffend())
  {
    std::cout << Name() << " clk: 0x" << std::hex << *iter1
              << " Gl1 clk: 0x" << *iter2 << std::dec << std::endl;
    iter1++;
    iter2++;
  }
  return;
}

int SingleTriggeredInput::ReadEvent()
{
  if (m_EventDeque.empty())
  {
    if (!EventAlignmentProblem())
    {
      std::cout << Name() << ":all events done" << std::endl;
      AllDone(1);
    }
    return -1;
  }
  if (Verbosity() > 1)
  {
    std::cout << "deque size: " << m_EventDeque.size() << std::endl;
  }
  m_FEMEventNrSet.clear();  // reset for new event, should be called at event cleanup
  Event *evt = m_EventDeque.front();
  m_EventDeque.pop_front();
  RunNumber(evt->getRunNumber());
  //  std::cout << "Handling event " << evt->getEvtSequence();
  //  CaloPacket *newhit = new CaloPacketv1();
  std::vector<Packet *> pktvec = evt->getPacketVector();
  std::vector<CaloPacket *> calopacketvector;
  for (auto *packet : pktvec)
  {
    int packet_id = packet->getIdentifier();
    CaloPacket *newhit = findNode::getClass<CaloPacket>(m_topNode, packet_id);
    calopacketvector.push_back(newhit);
    uint64_t packetbco = packet->lValue(0, "CLOCK");
    if (m_DitchPackets)
    {
      newhit->setStatus(OfflinePacket::PACKET_DROPPED);
      std::cout << Name() << " ditching packet " << packet_id << " from event " << evt->getEvtSequence() << std::endl;
      continue;
    }
    if (packet->getStatus())
    {
      newhit->setStatus(OfflinePacket::PACKET_CORRUPT);
      std::cout << Name() << " ditching corrupt packet " << packet_id << " from event " << evt->getEvtSequence() << std::endl;
      continue;
    }
    newhit->setStatus(OfflinePacket::PACKET_OK);
    newhit->setPacketEvtSequence(packet->iValue(0, "EVTNR"));
    int nr_modules = packet->iValue(0, "NRMODULES");
    int nr_channels = packet->iValue(0, "CHANNELS");
    int nr_samples = packet->iValue(0, "SAMPLES");
    newhit->setNrModules(nr_modules);
    newhit->setNrChannels(nr_channels);
    newhit->setNrSamples(nr_samples);
    newhit->setIdentifier(packet_id);
    newhit->setBCO(packetbco);
    //     std::cout << ", clock :" << packet->lValue(0, "CLOCK") << std::endl;
    for (int ifem = 0; ifem < nr_modules; ifem++)
    {
      newhit->setFemClock(ifem, packet->iValue(ifem, "FEMCLOCK"));
      newhit->setFemEvtSequence(ifem, packet->iValue(ifem, "FEMEVTNR"));
      newhit->setFemSlot(ifem, packet->iValue(ifem, "FEMSLOT"));
      newhit->setChecksumLsb(ifem, packet->iValue(ifem, "CHECKSUMLSB"));
      newhit->setChecksumMsb(ifem, packet->iValue(ifem, "CHECKSUMMSB"));
      newhit->setCalcChecksumLsb(ifem, packet->iValue(ifem, "CALCCHECKSUMLSB"));
      newhit->setCalcChecksumMsb(ifem, packet->iValue(ifem, "CALCCHECKSUMMSB"));
      newhit->setFemStatus(ifem, CaloPacket::FEM_OK);
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
    delete packet;
    int iret = FemEventNrClockCheck(newhit);
    if (iret < 0)
    {
      newhit->Reset();
    }
  }
  m_DitchPackets = false;
  if (m_FEMEventNrSet.size() > 1)
  {
    std::cout << "FEM event number mismatch among packets, aborting combining, from now on GL1 only" << std::endl;
    for (auto *iter : calopacketvector)
    {
      iter->Reset();
    }
    EventAlignmentProblem(1);
    for (auto *iter : m_EventDeque)
    {
      delete iter;
    }
    m_EventDeque.clear();
  }
  delete evt;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SingleTriggeredInput::checkfirstsebevent()
{
  // copy arrays into vectors for easier searching
  std::deque<uint64_t> gl1clkvector;
  for (const auto *iter = Gl1Input()->clkdiffbegin(); iter != Gl1Input()->clkdiffend(); ++iter)
  {
    gl1clkvector.push_back(*iter);
  }
  std::deque<uint64_t> myclkvector;
  for (const auto *iter = clkdiffbegin(); iter != clkdiffend(); ++iter)
  {
    myclkvector.push_back(*iter);
  }
  // remove matches from the front of the deques
  // the tricky part is that if there is a mismatch in clock diffs it is caused by
  // the previous channel, so we have to subtract the last one where the clock diff matched
  int badindex = -1;
  while (*(myclkvector.begin()) == *(gl1clkvector.begin()))
  {
    badindex++;
    //    std::cout << "removing " << std::hex << *(myclkvector.begin()) << std::dec << std::endl;
    myclkvector.pop_front();
    gl1clkvector.pop_front();
  }
  while (!myclkvector.empty())
  {
    auto it = std::search(gl1clkvector.begin(), gl1clkvector.end(), myclkvector.begin(), myclkvector.end());
    if (it == gl1clkvector.end())
    {
      if (Verbosity() > 0)
      {
        std::cout << "found no match comparing seb and Gl1 clock diffs" << std::endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "found matching sequence" << std::endl;
      }
      break;
    }
    myclkvector.pop_front();
  }
  if (myclkvector.empty())
  {
    std::cout << "cannot recover this problem - it is not a missing first event from the seb" << std::endl;
    EventAlignmentProblem(1);
  }
  // reshuffle clock refs
  for (unsigned int i = badindex; i < pooldepth - 1; i++)
  {
    m_bclkarray[i] = m_bclkarray[i + 1];
    m_bclkdiffarray[i] = m_bclkdiffarray[i + 1];
  }
  m_bclkarray[pooldepth] = std::numeric_limits<uint64_t>::max();
  //    std::cout << Name() << " deleting event " << m_EventDeque[badindex]->getEvtSequence() << std::endl;
  delete m_EventDeque[badindex];
  m_EventDeque.erase(m_EventDeque.begin() + badindex);
  // read next event to fill our arrays up to pooldepth
  Event *evt = GetEventIterator()->getNextEvent();
  while (!evt)
  {
    fileclose();
    if (!OpenNextFile())
    {
      FilesDone(1);
      if (Verbosity() > 0)
      {
        std::cout << "no more events to read, deque depth: " << m_EventDeque.size() << std::endl;
      }
      return -1;
    }
    evt = GetEventIterator()->getNextEvent();
  }
  m_EventsThisFile++;
  if (evt->getEvtType() != DATAEVENT)
  {
    std::cout << Name() << " things will crash and burn: non data event: " << evt->getEvtSequence() << std::endl;
    delete evt;
  }
  uint64_t myClock = GetClock(evt);
  m_bclkarray[pooldepth] = myClock;  // GetClock(evt);
  if (m_bclkarray[pooldepth] < m_bclkarray[pooldepth - 1])
  {
    m_bclkarray[pooldepth] += 0x100000000;
  }
  m_bclkdiffarray[pooldepth - 1] = m_bclkarray[pooldepth] - m_bclkarray[pooldepth - 1];
  m_Event++;
  m_EventDeque.push_back(evt);

  return 0;
}

void SingleTriggeredInput::RunCheck()
{
  bool isequal = std::equal(clkdiffbegin(), clkdiffend(), Gl1Input()->clkdiffbegin());
  if (!isequal)
  {
    std::cout << Name() << " and GL1 clock diffs differ (may be false positive at end of run), here is the dump:" << std::endl;
    dumpdeque();  // dump the clock diffs so we can see in the log
    const auto *iter1 = clkdiffbegin();
    const auto *iter2 = Gl1Input()->clkdiffbegin();
    if (*(++iter1) != *(++iter2) && firstclockcheck)
    {
      // this only works for the first event we process,
      // we know it is the first event we handle so setting
      // this flag is sufficient. Subsequent mismatches will show up as the second event mismatching
      // since we leave the last event in the deque
      std::cout << Name() << " first event problem, check subsequent events if our first event is off" << std::endl;
      iter1++;
      iter2++;
      if (std::equal(iter1, clkdiffend(), iter2))
      {
        std::cout << "subsequent events are good, first event is off reset packets from first event" << std::endl;
        m_DitchPackets = true;
      }
      else
      {
        std::cout << "check subsequent events failed, this is not the problem" << std::endl;
        if (checkfirstsebevent() == 0)
        {
          std::cout << Name() << " started with bad event" << std::endl;
        }
        else
        {
          std::cout << Name() << " first event in seb is not the problem either" << std::endl;
        }
      }
    }
    else
    {
      iter1 = clkdiffbegin();
      iter2 = Gl1Input()->clkdiffbegin();
      //        auto iter3 = beginclock();
      //        auto iter4 = Gl1Input()->beginclock();
      int position = 0;
      //        int ifirst = 1;
      while (iter1 != clkdiffend())
      {
        //           std::cout  << "bclk1: 0x" << *iter3 << ", gl1bclk: 0x" << *iter4 << std::dec << std::endl;
        // this catches the condition where there is no gl1 packet (14001)
        // the GL1 data sometimes has a corrupt last data event
        if (*iter1 != std::numeric_limits<uint64_t>::max())
        {
          std::cout << Name() << ": ";
          m_EventDeque[position]->identify();
        }
        if (*iter1 != *iter2)
        {
          if (*iter1 == std::numeric_limits<uint64_t>::max())
          {
            std::cout << Name() << " No more events found marked by uint64_t max clock" << std::endl;
            FilesDone(1);
            break;
          }
          if (Verbosity() > 0)
          {
            std::cout << "remove Event " << m_EventDeque[position]->getEvtSequence() << " clock: 0x"
                      << std::hex << GetClock(m_EventDeque[position]) << ", clkdiff(me) 0x"
                      << *iter1 << ", clkdiff(gl1) 0x" << *iter2
                      << std::dec << std::endl;
          }
          // position is the first bad index,
          for (size_t i = position; i < m_EventDeque.size(); ++i)
          {
            delete m_EventDeque[i];
          }
          m_EventDeque.erase(m_EventDeque.begin() + (position), m_EventDeque.end());
          break;
        }

        // std::cout <<  "good Event " << m_EventDeque[position]->getEvtSequence() << " clock: " << m_EventDeque[position]->getPacket(6067)->lValue(0, "CLOCK")<< std::endl ;

        ++position;
        ++iter1;
        ++iter2;
        // if (ifirst)
        // {
        //   ifirst = 0;
        // }
        // else
        // {
        //   ++iter3;
        //   ++iter4;
        // }
      }
      if (FilesDone() == 0)
      {
        std::cout << "Event Misalignment, processing remaining good events" << std::endl;
        EventAlignmentProblem(1);
      }
      //        FilesDone(1);
    }  // test for first event
  }  // end test for equality
  else
  {
    //	dumpdeque();
  }
  firstclockcheck = false;
  //	std::cout << "we are good" << std::endl;
  return;
}

bool SingleTriggeredInput::DoneFilling() const
{
  //  std::cout << Name() << " deq size: " << m_EventDeque.size() << std::endl;
  if (FilesDone() || m_EventDeque.size() >= pooldepth)
  {
    return true;
  }
  return false;
}

void SingleTriggeredInput::ResetClockDiffCounters()
{
  uint64_t tmp = m_bclkarray[pooldepth];
  m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkarray[0] = tmp;
  return;
}

uint64_t SingleTriggeredInput::calccurrentclockdiff(const int index, Event *evt)
{
  uint64_t myClock = GetClock(evt);
  uint64_t currentclockdiff{0};
  if (myClock < m_bclkarray[index])
  {
    currentclockdiff = myClock + 0x100000000 - m_bclkarray[index];
  }
  else
  {
    currentclockdiff = myClock - m_bclkarray[index];
  }
  if (Verbosity() > 1)
  {
    std::cout << Name() << " calc clkdiff: 0x" << std::hex << currentclockdiff
              << " gl1: 0x" << Gl1Input()->get_clkdiff(index) << std::dec << std::endl;
  }
  return currentclockdiff;
}

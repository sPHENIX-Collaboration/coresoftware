#include "SingleTriggerInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacket.h>
#include <phool/phool.h>

#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <set>
#include <utility>  // for pair
#include <vector>

SingleTriggerInput::SingleTriggerInput(const std::string &name)
  : Fun4AllBase(name)
{
}

SingleTriggerInput::~SingleTriggerInput()
{
  for (auto &openfiles : m_PacketDumpFile)
  {
    openfiles.second->close();
    delete openfiles.second;
  }
  m_PacketDumpFile.clear();
  delete m_EventIterator;
}

int SingleTriggerInput::fileopen(const std::string &filenam)
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

int SingleTriggerInput::fileclose()
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

void SingleTriggerInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto &bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "FEM: " << feeiter << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

bool SingleTriggerInput::CheckPoolDepth(const uint64_t bclk)
{
  // if (m_FEEBclkMap.size() < 10)
  // {
  //   std::cout << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
  //   return true;
  // }
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    if (iter.second < bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEE " << iter.first << " beamclock 0x" << std::hex << iter.second
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return true;
    }
  }
  return false;
}

void SingleTriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  m_BclkStack.erase(currentbclk);
  m_BeamClockFEE.erase(currentbclk);
  return;
}

void SingleTriggerInput::ddumppacket(Packet *pkt)
{
  int packetid = pkt->getIdentifier();
  if (m_PacketDumpFile.find(packetid) == m_PacketDumpFile.end())
  {
    std::string fname = "packet_" + std::to_string(packetid) + ".ddump";
    std::ofstream *dumpfile = new std::ofstream(fname);
    // dumpfile.open(fname);
    m_PacketDumpFile.insert(std::make_pair(packetid, dumpfile));
    m_PacketDumpCounter.insert(std::make_pair(packetid, m_ddump_flag));
  }
  if (m_PacketDumpCounter[packetid] != 0)
  {
    pkt->dump(*m_PacketDumpFile[packetid]);
    m_PacketDumpCounter[packetid]--;
  }
  return;
}

int SingleTriggerInput::EventNumberOffset(const int packetid)
{
  // initialize to zero, if map entry exists it will not overwrite it
  // just return false in retcode.second
  auto retcode = m_EventNumberOffset.insert(std::make_pair(packetid, m_DefaultEventNumberOffset));
  if (Verbosity() > 2)
  {
    if (retcode.second)
    {
      std::cout << PHWHERE << " Inserted " << m_DefaultEventNumberOffset << " as event offset for packet "
                << packetid << std::endl;
    }
  }
  return m_EventNumberOffset[packetid];
}

void SingleTriggerInput::AdjustEventNumberOffset(const int packetid, const int offset)
{
  if (m_EventNumberOffset.find(packetid) == m_EventNumberOffset.end())
  {
    return;
  }
  m_EventNumberOffset[packetid] += offset;
}

int SingleTriggerInput::AdjustPacketMap(int pktid, int evtoffset)
{
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << " adjusting local " << Name()
              << " packet map for packet " << pktid
              << " with offset " << evtoffset << std::endl;
  }
  std::vector<int> eventnumbers;
  for (auto packetmapiter = m_PacketMap.rbegin(); packetmapiter != m_PacketMap.rend(); ++packetmapiter)
  {
    eventnumbers.push_back(packetmapiter->first);
  }

  for (auto evtnumiter : eventnumbers)
  {
    int lastevent = evtnumiter;
    int newevent = lastevent + evtoffset;
    //    for (auto pktiter : m_PacketMap[lastevent])
    for (std::vector<OfflinePacket *>::iterator pktiter = m_PacketMap[lastevent].begin(); pktiter != m_PacketMap[lastevent].end(); ++pktiter)
    {
      if ((*pktiter)->getIdentifier() == pktid)
      {
        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << " need to move packet " << (*pktiter)->getIdentifier() << std::endl;
        }
        //      trivial variables give no speed benefit from using std::move
        //	m_PacketMap[newevent].push_back(std::move(*pktiter));
        m_PacketMap[newevent].push_back(*pktiter);
        m_PacketMap[lastevent].erase(pktiter);
        break;
      }
    }
  }
  return 0;
}

int SingleTriggerInput::AdjustEventOffset(int evtoffset)
{
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << " adjusting local " << Name()
              << " all packets with offset " << evtoffset << std::endl;
  }
  std::vector<int> eventnumbers;
  // needs separate cases so we don't overwrite existing entries
  if (evtoffset < 0)  // for negative offsets start at the beginning and move down (into empty space)
  {
    for (auto &packetmapiter : m_PacketMap)
    {
      eventnumbers.push_back(packetmapiter.first);
    }
  }
  else  // for positive offsets start at the end and move up (into empty space)
  {
    for (auto &packetmapiter : m_PacketMap)
    {
      eventnumbers.push_back(packetmapiter.first);
    }
  }

  for (auto evtnumiter : eventnumbers)
  {
    int lastevent = evtnumiter;
    int newevent = lastevent + evtoffset;
    //    for (auto pktiter : m_PacketMap[lastevent])
    for (std::vector<OfflinePacket *>::iterator pktiter = m_PacketMap[lastevent].begin(); pktiter != m_PacketMap[lastevent].end(); ++pktiter)
    {
      //      if ((*pktiter)->getIdentifier() == pktid)
      {
        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << " need to move packet " << (*pktiter)->getIdentifier() << std::endl;
        }
        //      trivial variables give no speed benefit from using std::move
        //	m_PacketMap[newevent].push_back(std::move(*pktiter));
        m_PacketMap[newevent].push_back(*pktiter);
        m_PacketMap[lastevent].erase(pktiter);
        break;
      }
    }
  }
  for (auto evtnumiter : m_EventNumberOffset)
  {
    evtnumiter.second += evtoffset;
  }
  return 0;
}

bool SingleTriggerInput::GetSomeMoreEvents(const unsigned int keep)
{
  if (AllDone())
  {
    return false;
  }
  if (m_PacketMap.empty())
  {
    return true;
  }
  if (Verbosity() > 21)
  {
    std::cout << PHWHERE << Name() << ": first event: " << m_PacketMap.begin()->first
              << " last event: " << m_PacketMap.rbegin()->first << " size: " << m_PacketMap.size()
              << ", keep: " << keep
              << std::endl;
  }
// how many events should be stored upstream (keep) plus number of events kept locally
  if (m_PacketMap.size() < std::max(2U, keep + m_LocalPoolDepth))  // at least 2 events in pool
  {
    return true;
  }
  return false;
}

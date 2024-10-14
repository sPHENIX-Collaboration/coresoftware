#include "SingleTriggerInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacket.h>
#include <phool/phool.h>

#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <TSystem.h>

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

int SingleTriggerInput::SetFEMEventRefPacketId(const int pktid)
{
  if (m_FEMEventRefPacketId > 0)
  {
    return m_FEMEventRefPacketId;
  }
// from https://indico.bnl.gov/event/24287/contributions/94631/attachments/56257/96280/Software_20240730.pdf
// /* 6067 copying to 6068 6071 6072 6075 6076 6079 6080 */ seb00
// /* 6083 copying to 6084 6087 6088 6091 6092 6095 6096 */ seb01
// /* 6115 copying to 6116 6119 6120 6123 6124 6127 6128 */ seb02
// /* 6099 copying to 6100 6103 6104 6107 6108 6111 6112 */ seb03
// /* 6035 copying to 6036 6039 6040 6043 6044 6047 6048 */ seb04
// /* 6051 copying to 6052 6055 6056 6059 6060 6063 6064 */ seb05
// /* 6019 copying to 6020 6023 6024 6027 6028 6031 6032 */ seb06
// /* 6003 copying to 6004 6007 6008 6011 6012 6015 6016 */ seb07
// /* 6065 copying to 6066 6069 6070 6073 6074 6077 6078 */ seb08
// /* 6081 copying to 6082 6085 6086 6089 6090 6093 6094 */ seb09
// /* 6121 copying to 6122 6125 6126 6113 6114 6117 6118 */ seb10
// /* 6105 copying to 6106 6109 6110 6097 6098 6101 6102 */ seb11
// /* 6033 copying to 6034 6037 6038 6041 6042 6045 6046 */ seb12
// /* 6049 copying to 6050 6053 6054 6057 6058 6061 6062 */ seb13
// /* 6025 copying to 6026 6029 6030 6017 6018 6021 6022 */ seb14
// /* 6009 copying to 6010 6013 6014 6001 6002 6005 6006 */ seb15
// /* 8001 copying to 8002 8007 8008 7001 7002 7007 7008  */ seb16
// /* 8003 copying to 8004 8005 8006 7003 7004 7005 7006 */ seb17
// /* 9003 copying to 9002 9001 9006 9005 9004 */ seb20 -- sepd only
  static std::map<int, int> refpacketmap = {
    {6067, 6067}, {6068, 6067}, {6071, 6067}, {6072, 6067}, {6075, 6067}, {6076, 6067}, {6079, 6067}, {6080, 6067},
    {6083, 6083}, {6084, 6083}, {6087, 6083}, {6088, 6083}, {6091, 6083}, {6092, 6083}, {6095, 6083}, {6096, 6083},
    {6115, 6115}, {6116, 6115}, {6119, 6115}, {6120, 6115}, {6123, 6115}, {6124, 6115}, {6127, 6115}, {6128, 6115},
    {6099, 6099}, {6100, 6099}, {6103, 6099}, {6104, 6099}, {6107, 6099}, {6108, 6099}, {6111, 6099}, {6112, 6099},
    {6035, 6035}, {6036, 6035}, {6039, 6035}, {6040, 6035}, {6043, 6035}, {6044, 6035}, {6047, 6035}, {6048, 6035},
    {6051, 6051}, {6052, 6051}, {6055, 6051}, {6056, 6051}, {6059, 6051}, {6060, 6051}, {6063, 6051}, {6064, 6051},
    {6019, 6019}, {6020, 6019}, {6023, 6019}, {6024, 6019}, {6027, 6019}, {6028, 6019}, {6031, 6019}, {6032, 6019},
    {6003, 6003}, {6004, 6003}, {6007, 6003}, {6008, 6003}, {6011, 6003}, {6012, 6003}, {6015, 6003}, {6016, 6003},
    {6065, 6065}, {6066, 6065}, {6069, 6065}, {6070, 6065}, {6073, 6065}, {6074, 6065}, {6077, 6065}, {6078, 6065},
    {6081, 6081}, {6082, 6081}, {6085, 6081}, {6086, 6081}, {6089, 6081}, {6090, 6081}, {6093, 6081}, {6094, 6081},
    {6121, 6121}, {6122, 6121}, {6125, 6121}, {6126, 6121}, {6113, 6121}, {6114, 6121}, {6117, 6121}, {6118, 6121},
    {6105, 6105}, {6106, 6105}, {6109, 6105}, {6110, 6105}, {6097, 6105}, {6098, 6105}, {6101, 6105}, {6102, 6105},
    {6033, 6033}, {6034, 6033}, {6037, 6033}, {6038, 6033}, {6041, 6033}, {6042, 6033}, {6045, 6033}, {6046, 6033},
    {6049, 6049}, {6050, 6049}, {6053, 6049}, {6054, 6049}, {6054, 6049}, {6058, 6049}, {6061, 6049}, {6062, 6049},
    {6025, 6025}, {6026, 6025}, {6029, 6025}, {6030, 6025}, {6017, 6025}, {6018, 6025}, {6021, 6025}, {6022, 6025},
    {6009, 6009}, {6010, 6009}, {6013, 6009}, {6014, 6009}, {6001, 6009}, {6002, 6009}, {6005, 6009}, {6006, 6009},
    {8001, 8001}, {8002, 8001}, {8007, 8001}, {8008, 8001}, {7001, 8001}, {7002, 8001}, {7007, 8001}, {7008, 8001},
    {8003, 8003}, {8004, 8003}, {8005, 8003}, {7003, 8003}, {7004, 8003}, {7005, 8003}, {7006, 8003},
    {9003, 9003}, {9002, 9003}, {9001, 9003}, {9006, 9003}, {9004, 9003}
};
  auto refpacket_pair = refpacketmap.find(pktid);
  if (refpacket_pair != refpacketmap.end())
  {
    m_FEMEventRefPacketId = refpacket_pair->second;
    return m_FEMEventRefPacketId;
  }
  std::cout << PHWHERE << " could not find " << pktid << " in refpacket map" << std::endl;
  gSystem->Exit(1);
  exit(1);
}

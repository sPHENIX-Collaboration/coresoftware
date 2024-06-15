#include "SingleTriggerInput.h"

#include <frog/FROG.h>

#include <phool/phool.h>

#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <set>
#include <utility>  // for pair

SingleTriggerInput::SingleTriggerInput(const std::string &name)
  : Fun4AllBase(name)
{
}

SingleTriggerInput::~SingleTriggerInput()
{
  for (auto &openfiles : m_PacketDumpFile)
  {
    openfiles.second->close();
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
    //dumpfile.open(fname);
      m_PacketDumpFile.insert(std::make_pair(packetid,dumpfile));
      m_PacketDumpCounter.insert(std::make_pair(packetid,m_ddump_flag));
  }
  if (m_PacketDumpCounter[packetid] != 0)
  {
    pkt->dump(*m_PacketDumpFile[packetid]);
    m_PacketDumpCounter[packetid]--;
  }
  return;
}

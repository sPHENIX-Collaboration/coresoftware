#include "SingleEvtInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <frog/FROG.h>

#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <set>

SingleEvtInput::SingleEvtInput(const std::string &name, Fun4AllEvtInputPoolManager *inman)
  : Fun4AllBase(name)
  , m_InputMgr(inman)
{
  plist = new Packet *[100];
}

SingleEvtInput::~SingleEvtInput()
{
  delete m_EventIterator;
  delete[] plist;
}

void SingleEvtInput::FillPool(const unsigned int nbclks)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (m_EventIterator == nullptr)  // at startup this is a null pointer
  {
    OpenNextFile();
  }
  std::set<uint64_t> saved_beamclocks;
  do
  {
    Event *evt = m_EventIterator->getNextEvent();
    std::cout << "Fetching next Event" << std::endl;
    if (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt = m_EventIterator->getNextEvent();
      if (!evt)
      {
        std::cout << PHWHERE << "Event is nullptr" << std::endl;
        AllDone(1);
        return;
      }
    }
    m_RunNumber = evt->getRunNumber();
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 100);

    if (npackets == 100)
    {
      exit(1);
    }
    for (int i = 0; i < npackets; i++)
    {
      plist[i]->identify();
      //        if (plist[i]->iValue(0, "EVENCHECKSUMOK") != 0 && plist[i]->iValue(0, "ODDCHECKSUMOK") != 0)
      {
        int num_hits = plist[i]->iValue(0, "NR_HITS");
	std::set<uint64_t> bclk_set;
	for (int j=0; j<num_hits; j++)
	{
	  int FEE =  plist[i]->iValue(j, "FEE");
	  uint64_t gtm_bco = plist[i]->lValue(j, "BCO") + m_Rollover[FEE];
	  bclk_set.insert(gtm_bco);
	  if (gtm_bco < m_PreviousClock[FEE])
	  {
	    m_Rollover[FEE] += 0x100000000;
	    gtm_bco += m_Rollover[FEE];  // rollover makes sure our bclks are ascending even if we roll over the 40 bit counter
	  }
	  m_PreviousClock[FEE] = gtm_bco;
	  m_BeamClockFEE[gtm_bco].insert(FEE);
          m_FEEBclkMap[FEE] = gtm_bco;
	  std::cout << "evtno: " << EventSequence
		    << ", nr_hits: " << num_hits
		    << ", bco: 0x" << std::hex << gtm_bco << std::dec 
		    << ", FEE: " << FEE << std::endl;
	  // dummy check for the first event which is the problem for the calorimeters
	  // it is the last event from the previous run, so it's event number is > 0
	  // if (evtno > EventSequence)
	  // {
	  //   delete plist[i];
	  //   plist[i] = nullptr;
	  //   continue;
	  // }
	  plist[i]->convert();
	  // calculate "real" event number
	  // special events are counted, so the packet event counter is never the
	  // Event Sequence (bc the begin run event)
	  // also our packets are just 16bit counters, so we need to add the upper bits
	  // from the event sequence
	  // and our packet counters start at 0, while our events start at 1

	  //          evtno += m_EventNumberOffset + m_NumSpecialEvents + (EventSequence & 0xFFFF0000);
	}
	uint64_t latestbclk = 0;
	if (bclk_set.empty())
	{
	  delete plist[i];
	}
	else
	{
	for (auto bco_iter: bclk_set)
	{
	  saved_beamclocks.insert(bco_iter);
	  if (m_InputMgr)
	  {
	    m_InputMgr->AddPacket(bco_iter, plist[i]);
	  }
          else
	  {
	    m_BclkStack.insert(bco_iter);
	  }
	}
        latestbclk = *bclk_set.rbegin();
        auto packetvector = m_PacketStorageMap.find(latestbclk);
	if (packetvector == m_PacketStorageMap.end())
	{
	  std::vector<Packet *> packets;
	m_PacketStorageMap[latestbclk] = packets;
        packetvector =  m_PacketStorageMap.find(latestbclk);
	}
	packetvector->second.push_back(plist[i]);
	}
      }
      // else
      // {
      //   delete plist[i];
      // }
    }
    delete evt;
  } while (CheckPoolDepth(m_PacketStorageMap.begin()->first));
}

int SingleEvtInput::fileopen(const std::string &filenam)
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

int SingleEvtInput::fileclose()
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

void SingleEvtInput::Print(const std::string &what)
{
  if (what == "ALL" || what == "FEE")
  {
    for (auto bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter:  bcliter.second)
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
  if (what == "ALL" || what == "STORAGE")
  {
    for (auto bcliter : m_PacketStorageMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter:  bcliter.second)
      {
	std::cout << "Packet: " << feeiter->getIdentifier()
		  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
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

void SingleEvtInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (auto iter : m_PacketStorageMap)
  {
    if (iter.first <= bclk)
    {
      for (auto pktiter: iter.second)
      {
	delete pktiter;
      }
      toclearbclk.push_back(iter.first);
    }
    else
    {
      break;
    }
  }
  for (auto iter : toclearbclk)
  {
      m_PacketStorageMap.erase(iter);
  }
}
bool SingleEvtInput::CheckPoolDepth(const uint64_t bclk)
{
  if (m_FEEBclkMap.size() < 14)
  {
    std::cout << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
    return true;
  }
  for (auto iter : m_FEEBclkMap)
  {
    std::cout << "my bclk 0x" << std::hex << iter.second
	      << " req: 0x" << bclk << std::dec << std::endl;
    if (iter.second < bclk)
    {
      std::cout << "my beamclock 0x" << std::hex << iter.second
		<< " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      return true;
    }
  }
  return false;
}

void SingleEvtInput::ClearCurrentEvent()
{
// called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  m_BclkStack.erase(currentbclk);
  m_BeamClockFEE.erase(currentbclk);
  return;
}

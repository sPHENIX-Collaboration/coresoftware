#include "SinglePrdfInput.h"

#include "Fun4AllPrdfInputPoolManager.h"

#include <frog/FROG.h>

#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

SinglePrdfInput::SinglePrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman)
  : Fun4AllBase(name)
  , m_InputMgr(inman)
{
  plist = new Packet *[100];
  m_PacketEventNumberOffset = new int[100]{};
  //  std::fill_n(m_PacketEventNumberOffset, 100, 0);
}

SinglePrdfInput::~SinglePrdfInput()
{
  delete m_EventIterator;
  delete[] plist;
  delete[] m_PacketEventNumberOffset;
}

void SinglePrdfInput::FillPool(const unsigned int nevents)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (m_EventIterator == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }
  for (unsigned int ievt = 0; ievt < nevents; ievt++)
  {
    Event *evt = m_EventIterator->getNextEvent();
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
    if (Verbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      delete evt;
      continue;  // need handling for non data events
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 100);
    if (npackets == 100)
    {
      exit(1);
    }
    for (int i = 0; i < npackets; i++)
    {
      if (plist[i]->iValue(0, "EVENCHECKSUMOK") != 0 && plist[i]->iValue(0, "ODDCHECKSUMOK") != 0)
      {
        int evtno = plist[i]->iValue(0, "EVTNR");
        unsigned int bclk = plist[i]->iValue(0, "CLOCK");
        if (Verbosity() > 1)
        {
          std::cout << "packet " << plist[i]->getIdentifier() << " evt: " << evtno
                    << std::hex << " clock: 0x" << bclk << std::dec << std::endl;
        }
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
        evtno += m_EventNumberOffset + m_PacketEventNumberOffset[i] + m_NumSpecialEvents + (EventSequence & 0xFFFF0000);
        m_PacketMap[bclk].push_back(plist[i]);
        m_EvtSet.insert(evtno);
        m_Event.emplace_back(std::make_pair(evtno, bclk));
      }
      else
      {
        delete plist[i];
      }
    }
    // here we have all packets of a given event in our maps/vectors
    // first pass - check if beam clocks are identical
    if (Verbosity() > 1)
    {
      std::cout << "pktmap size : " << m_PacketMap.size() << std::endl;
      std::cout << "evt set size : " << m_EvtSet.size() << std::endl;
    }
    int common_event_number = *(m_EvtSet.begin());
    int common_beam_clock = m_PacketMap.begin()->first;
    if (m_PacketMap.size() == 1)  // all packets from the same beam clock
    {
      if (m_EvtSet.size() == 1)
      {
        if (Verbosity() > 1)
        {
          std::cout << "we are good evtno: " << *(m_EvtSet.begin())
                    << ", clock: " << m_PacketMap.begin()->first << std::endl;
        }
      }
      else
      {
        if (Verbosity() > 1)
        {
          std::cout << "We have multiple event numbers for bclk: 0x" << std::hex
                    << m_PacketMap.begin()->first << std::dec << std::endl;
          for (auto iter : m_EvtSet)
          {
            std::cout << "Event " << iter << std::endl;
          }
        }
        common_event_number = majority_eventnumber();
        if (Verbosity() > 1)
        {
          std::cout << "picked event no " << common_event_number << std::endl;
        }
        adjust_eventnumber_offset(common_event_number);
      }
      for (auto const &iter : m_PacketMap)
      {
        for (auto const &pktiter : iter.second)
        {
          m_InputMgr->AddPacket(common_event_number, pktiter);
        }
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "We have multiple beam clocks per event" << std::endl;
      }
      if (m_EvtSet.size() == 1)
      {
        if (Verbosity() > 1)
        {
          std::cout << "we are good evtno: " << *(m_EvtSet.begin())
                    << ", clock: " << m_PacketMap.begin()->first << std::endl;
        }
      }
      else
      {
        if (Verbosity() > 1)
        {
          std::cout << "We have multiple event numbers for bclk: 0x" << std::hex
                    << m_PacketMap.begin()->first << std::dec << std::endl;
          for (auto iter : m_EvtSet)
          {
            std::cout << "Event " << iter << std::endl;
          }
        }
        common_event_number = majority_eventnumber();
        if (Verbosity() > 1)
        {
          std::cout << "picked event no " << common_event_number << std::endl;
        }
        adjust_eventnumber_offset(common_event_number);
      }
      common_beam_clock = majority_beamclock();
      if (Verbosity() > 1)
      {
        std::cout << "picked bclk: " << std::hex << common_beam_clock << std::dec << std::endl;
      }
      // for time being clean out packets which do not match
      for (auto const &iter : m_PacketMap)
      {
        for (auto pktiter : iter.second)
        {
          if (pktiter->iValue(0, "CLOCK") == common_beam_clock)
          {
            if (Verbosity() > 1)
            {
              std::cout << "adding packet " << pktiter->getIdentifier() << " beam clock "
                        << std::hex << pktiter->iValue(0, "CLOCK") << std::dec << std::endl;
            }
            m_InputMgr->AddPacket(common_event_number, pktiter);
          }
          else
          {
            if (Verbosity() > 1)
            {
              std::cout << "Deleting packet " << pktiter->getIdentifier() << " beam clock "
                        << std::hex << pktiter->iValue(0, "CLOCK") << " common bclk: "
                        << common_beam_clock << std::dec << std::endl;
            }
            m_InputMgr->UpdateDroppedPacket(pktiter->getIdentifier());
            delete pktiter;
          }
        }
      }
    }
    m_InputMgr->AddBeamClock(common_event_number, common_beam_clock, this);
    if (m_MeReferenceFlag)
    {
      m_InputMgr->SetReferenceClock(common_event_number, common_beam_clock);
    }
    m_PacketMap.clear();
    m_EvtSet.clear();
    m_Event.clear();
    delete evt;
  }
}

void SinglePrdfInput::adjust_eventnumber_offset(const int decided_evtno)
{
  for (unsigned int i = 0; i < m_Event.size(); i++)
  {
    if (m_Event[i].first != decided_evtno)
    {
      m_PacketEventNumberOffset[i] -= (m_Event[i].first - decided_evtno);
      if (Verbosity() > 1)
      {
        std::cout << "my evtno: " << m_Event[i].first << ", decided: " << decided_evtno
                  << ", adjustment: " << m_Event[i].first - decided_evtno << std::endl;
        std::cout << "adjusting event number offset for " << i << " to " << m_PacketEventNumberOffset[i] << std::endl;
      }
    }
  }
}

int SinglePrdfInput::majority_eventnumber()
{
  std::map<int, int> evtcnt;
  for (auto iter : m_Event)
  {
    evtcnt[iter.first]++;
  }
  int imax = -1;
  int evtno = -1;
  for (auto iter : evtcnt)
  {
    if (iter.second > imax)
    {
      evtno = iter.first;
      imax = iter.second;
    }
  }
  return evtno;
}

int SinglePrdfInput::majority_beamclock()
{
  std::map<int, int> evtcnt;
  for (auto iter : m_Event)
  {
    evtcnt[iter.second]++;
    if (Verbosity() > 1)
    {
      std::cout << "adding clk: " << std::hex << iter.second << std::dec
                << " current counter: " << evtcnt[iter.second] << std::endl;
    }
  }
  int imax = -1;
  int bclk = INT_MAX;
  for (auto iter : evtcnt)
  {
    if (iter.second > imax)
    {
      bclk = iter.first;
      imax = iter.second;
    }
  }
  return bclk;
}

int SinglePrdfInput::fileopen(const std::string &filenam)
{
  std::cout << PHWHERE << Name() << ": trying to open " << filenam << std::endl;
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

int SinglePrdfInput::fileclose()
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

void SinglePrdfInput::MakeReference(const bool b)
{
  m_MeReferenceFlag = b;
  if (b && m_InputMgr)
  {
    m_InputMgr->SetReferenceInputMgr(this);
  }
  return;
}

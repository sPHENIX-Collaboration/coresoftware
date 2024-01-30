#include "SingleGl1PrdfInput.h"

#include "Fun4AllPrdfInputPoolManager.h"
#include "Fun4AllPrdfInputTriggerManager.h"

#include <frog/FROG.h>

#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

SingleGl1PrdfInput::SingleGl1PrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman)
  : SinglePrdfInput(name, inman)
{
  plist = new Packet *[100];
  m_PacketEventNumberOffset = new int[100]{};
}

SingleGl1PrdfInput::SingleGl1PrdfInput(const std::string &name, Fun4AllPrdfInputTriggerManager *inman)
  : SinglePrdfInput(name, inman)
{
  plist = new Packet *[100];
  m_PacketEventNumberOffset = new int[100]{};
}

SingleGl1PrdfInput::~SingleGl1PrdfInput()
{
  delete[] plist;
  delete[] m_PacketEventNumberOffset;
}

void SingleGl1PrdfInput::FillPool(const unsigned int nevents)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventIterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }
  for (unsigned int ievt = 0; ievt < nevents; ievt++)
  {
    Event *evt = GetEventIterator()->getNextEvent();
    if (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt = GetEventIterator()->getNextEvent();
      if (!evt)
      {
        std::cout << PHWHERE << "Event is nullptr" << std::endl;
        AllDone(1);
        return;
      }
    }
    RunNumber(evt->getRunNumber());
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
    //    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 100);
    if (npackets == 100)
    {
      exit(1);
    }
    int evtno = evt->getEvtSequence();
    for (int i = 0; i < npackets; i++)
    {
      //     if (plist[i]->iValue(0, "CHECKSUMOK") != 0)
      if (plist[i])
      {
        //        int evtno = plist[i]->iValue(0, "EVTNR");
        uint64_t bclk = plist[i]->lValue(0, "BCO");
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
        // calculate "real" event number
        // special events are counted, so the packet event counter is never the
        // Event Sequence (bc the begin run event)
        // also our packets are just 16bit counters, so we need to add the upper bits
        // from the event sequence
        // and our packet counters start at 0, while our events start at 1
        //        evtno += m_EventNumberOffset + m_PacketEventNumberOffset[i] + m_NumSpecialEvents;
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
    uint64_t common_beam_clock = m_PacketMap.begin()->first;
    if (m_PacketMap.size() == 1)  // all packets from the same beam clock
    {
      if (m_EvtSet.size() == 1)
      {
        if (Verbosity() > 1)
        {
          std::cout << "we are good evtno: " << *(m_EvtSet.begin())
                    << std::hex << ", clock: 0x" << m_PacketMap.begin()->first
                    << std::dec << std::endl;
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
          if (InputMgr())
          {
            InputMgr()->AddPacket(common_event_number, pktiter);
          }
          else
          {
            TriggerInputMgr()->AddPacket(common_event_number, pktiter);
          }
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
          if (((uint64_t) pktiter->lValue(0, "BCO")) == common_beam_clock)
          {
            if (Verbosity() > 1)
            {
              std::cout << "adding packet " << pktiter->getIdentifier() << " beam clock "
                        << std::hex << pktiter->lValue(0, "BCO") << std::dec << std::endl;
            }
            if (InputMgr())
            {
              InputMgr()->AddPacket(common_event_number, pktiter);
            }
            else
            {
              TriggerInputMgr()->AddPacket(common_event_number, pktiter);
            }
          }
          else
          {
            if (Verbosity() > 1)
            {
              std::cout << "Deleting packet " << pktiter->getIdentifier() << " beam clock "
                        << std::hex << pktiter->iValue(0, "CLOCK") << " common bclk: "
                        << common_beam_clock << std::dec << std::endl;
            }
            if (InputMgr())
            {
              InputMgr()->UpdateDroppedPacket(pktiter->getIdentifier());
            }
            else
            {
              TriggerInputMgr()->UpdateDroppedPacket(pktiter->getIdentifier());
            }

            delete pktiter;
          }
        }
      }
    }
    if (InputMgr())
    {
      InputMgr()->AddBeamClock(common_event_number, common_beam_clock, this);
    }
    else
    {
      TriggerInputMgr()->AddBeamClock(common_event_number, common_beam_clock, this);
    }

    if (ReferenceFlag())
    {
      if (InputMgr())
      {
        InputMgr()->SetReferenceClock(common_event_number, common_beam_clock);
      }
      else
      {
        TriggerInputMgr()->SetReferenceClock(common_event_number, common_beam_clock);
      }
    }
    m_PacketMap.clear();
    m_EvtSet.clear();
    m_Event.clear();
    delete evt;
  }
}

void SingleGl1PrdfInput::adjust_eventnumber_offset(const int decided_evtno)
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

int SingleGl1PrdfInput::majority_eventnumber()
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

int SingleGl1PrdfInput::majority_beamclock()
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

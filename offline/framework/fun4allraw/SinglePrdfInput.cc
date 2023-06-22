#include "SinglePrdfInput.h"

#include "Fun4AllPrdfInputPoolManager.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

SinglePrdfInput::SinglePrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman)
  : Fun4AllBase(name)
  , m_InputMgr(inman)
{
  plist = new Packet *[100];
}

SinglePrdfInput::~SinglePrdfInput()
{
  delete m_EventIterator;
}

void SinglePrdfInput::AddPrdfInputFile(const std::string &filename)
{
  m_FileList.push_back(filename);
  m_FileListCopy.push_back(filename);
}

void SinglePrdfInput::FillPool()
{
  if (m_EventMap.size() > m_LowWaterMark)
  {
    return;
  }

  if (m_EventIterator == nullptr)
  {
    int status = 0;
    std::cout << "opening " << m_FileList.front() << std::endl;
    m_EventIterator = new fileEventiterator(m_FileList.front().c_str(), status);
    if (status)
    {
      delete m_EventIterator;
      m_EventIterator = nullptr;
      return;
    }
  }
  if (m_PoolEvents <= m_LowWaterMark)
  {
    while (m_PoolEvents <= m_PoolDepth)
    {
      Event *evt = m_EventIterator->getNextEvent();
      m_RunNumber = evt->getRunNumber();
      evt->identify();
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
        if (plist[i]->iValue(0, "EVENCHECKSUMOK") != 0 && plist[i]->iValue(0, "ODDCHECKSUMOK") != 0)
        {
          int evtno = plist[i]->iValue(0, "EVTNR");
          // dummy check for the first event which is the problem for the calorimeters
          // it is the last event from the previous run, so it's event number is > 0
          if (evtno > EventSequence)
          {
            delete plist[i];
            plist[i] = nullptr;
            continue;
          }
          plist[i]->convert();
          // calculate "real" event number
          // special events are counted, so the packet event counter is never the
          // Event Sequence (bc the begin run event)
          // also our packets are just 16bit counters, so we need to add the upper bits
          // from the event sequence
          // and our packet counters start at 0, while our events start at 1

          evtno += m_EventNumberOffset + m_NumSpecialEvents + (EventSequence & 0xFFFF0000);
          m_InputMgr->AddPacket(evtno, plist[i]);
        }
      }
      m_PoolEvents++;
    }
  }
}

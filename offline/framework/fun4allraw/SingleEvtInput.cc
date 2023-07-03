#include "SingleEvtInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <frog/FROG.h>

#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

SingleEvtInput::SingleEvtInput(const std::string &name, Fun4AllEvtInputPoolManager *inman)
  : Fun4AllBase(name)
  , m_InputMgr(inman)
{
  plist = new Packet *[100];
}

SingleEvtInput::~SingleEvtInput()
{
  delete m_EventIterator;
  delete [] plist;
}

void SingleEvtInput::FillPool(const unsigned int nevents)
{
  if (AllDone()) // no more files and all events read
  {
    return;
  }
  while (m_EventIterator == nullptr) // at startup this is a null pointer
  {
    OpenNextFile();
  }
  for (unsigned int ievt=0; ievt<nevents; ievt++)
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
	// else
	// {
	//   delete plist[i];
	// }
      }
      delete evt;
  }
  
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


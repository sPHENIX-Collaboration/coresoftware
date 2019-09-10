#include "Fun4AllPrdfInputManager.h"

#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "Fun4AllSyncManager.h"
#include "Fun4AllUtils.h"

#include <ffaobjects/SyncObject.h>    // for SyncObject
#include <ffaobjects/SyncObjectv1.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>             // for PHNode
#include <phool/PHNodeIterator.h>     // for PHNodeIterator
#include <phool/phool.h>              // for PHWHERE

#include <Event/Event.h>
#include <Event/Eventiterator.h>      // for Eventiterator
#include <Event/fileEventiterator.h>

#include <cassert>
#include <cstdlib>
#include <iostream>                   // for operator<<, basic_ostream, endl
#include <utility>                    // for pair

using namespace std;

Fun4AllPrdfInputManager::Fun4AllPrdfInputManager(const string &name, const string &prdfnodename, const string &topnodename)
  : Fun4AllInputManager(name, prdfnodename, topnodename)
  , m_Segment(-999)
  , m_EventsTotal(0)
  , m_EventsThisFile(0)
  , m_Event(nullptr)
  , m_SaveEvent(nullptr)
  , m_EventIterator(nullptr)
  , m_SyncObject(new SyncObjectv1())
  , m_PrdfNodeName(prdfnodename)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());
  PHNodeIterator iter(m_topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  if (!PrdfNode)
  {
    PHDataNode<Event> *newNode = new PHDataNode<Event>(m_Event, m_PrdfNodeName, "Event");
    m_topNode->addNode(newNode);
  }
  return;
}

Fun4AllPrdfInputManager::~Fun4AllPrdfInputManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete m_SyncObject;
}

int Fun4AllPrdfInputManager::fileopen(const string &filenam)
{
  if (IsOpen())
  {
    cout << "Closing currently open file "
         << FileName()
         << " and opening " << filenam << endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  string fname = frog.location(FileName());
  if (Verbosity() > 0)
  {
    cout << Name() << ": opening file " << FileName() << endl;
  }
  int status = 0;
  m_EventIterator = new fileEventiterator(fname.c_str(), status);
  m_EventsThisFile = 0;
  if (status)
  {
    delete m_EventIterator;
    m_EventIterator = nullptr;
    cout << PHWHERE << Name() << ": could not open file " << fname << endl;
    return -1;
  }
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(fname);
  m_Segment = runseg.second;
  IsOpen(1);
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int Fun4AllPrdfInputManager::run(const int nevents)
{
readagain:
  if (!IsOpen())
  {
    if (FileListEmpty())

    {
      if (Verbosity() > 0)
      {
        cout << Name() << ": No Input file open" << endl;
      }
      return -1;
    }
    else
    {
      if (OpenNextFile())
      {
        cout << Name() << ": No Input file from filelist opened" << endl;
        return -1;
      }
    }
  }
  if (Verbosity() > 3)
  {
    cout << "Getting Event from " << Name() << endl;
  }
  //  cout << "running event " << nevents << endl;
  PHNodeIterator iter(m_topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  if (m_SaveEvent)  // if an event was pushed back, copy saved pointer and reset m_SaveEvent pointer
  {
    m_Event = m_SaveEvent;
    m_SaveEvent = nullptr;
    m_EventsThisFile--;
    m_EventsTotal--;
  }
  else
  {
    m_Event = m_EventIterator->getNextEvent();
  }
  PrdfNode->setData(m_Event);
  if (!m_Event)
  {
    fileclose();
    goto readagain;
  }
  if (Verbosity() > 1)
  {
    cout << Name() << " PRDF run " << m_Event->getRunNumber() << ", evt no: " << m_Event->getEvtSequence() << endl;
  }
  m_EventsTotal++;
  m_EventsThisFile++;
  SetRunNumber(m_Event->getRunNumber());
  MySyncManager()->PrdfEvents(m_EventsThisFile);
  MySyncManager()->SegmentNumber(m_Segment);
  MySyncManager()->CurrentEvent(m_Event->getEvtSequence());
  m_SyncObject->EventCounter(m_EventsThisFile);
  m_SyncObject->SegmentNumber(m_Segment);
  m_SyncObject->RunNumber(m_Event->getRunNumber());
  m_SyncObject->EventNumber(m_Event->getEvtSequence());
  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
  {
    ResetEvent();
    goto readagain;
  }
  return 0;
}

int Fun4AllPrdfInputManager::fileclose()
{
  if (!IsOpen())
  {
    cout << Name() << ": fileclose: No Input file open" << endl;
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

void Fun4AllPrdfInputManager::Print(const string &what) const
{
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllPrdfInputManager::ResetEvent()
{
  PHNodeIterator iter(m_topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  PrdfNode->setData(nullptr);  // set pointer in Node to nullptr before deleting it
  delete m_Event;
  m_Event = nullptr;
  m_SyncObject->Reset();
  return 0;
}

int Fun4AllPrdfInputManager::PushBackEvents(const int i)
{
  // PushBackEvents is supposedly pushing events back on the stack which works
  // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
  // A special case is when the synchronization fails and we need to only push back a single
  // event. In this case we save the m_Event pointer as m_SaveEvent which is used in the run method
  // instead of getting the next event.
  if (i > 0)
  {
    if (i == 1 && m_Event)  // check on m_Event pointer makes sure it is not done from the cmd line
    {
      m_SaveEvent = m_Event;
      return 0;
    }
    cout << PHWHERE << Name()
         << " Fun4AllPrdfInputManager cannot push back " << i << " events into file"
         << endl;
    return -1;
  }
  if (!m_EventIterator)
  {
    cout << PHWHERE << Name()
         << " no file open" << endl;
    return -1;
  }
  // Skipping events is implemented as
  // pushing a negative number of events on the stack, so in order to implement
  // the skipping of events we read -i events.
  int nevents = -i;  // negative number of events to push back -> skip num events
  int errorflag = 0;
  while (nevents > 0 && !errorflag)
  {
    m_Event = m_EventIterator->getNextEvent();
    if (!m_Event)
    {
      cout << "Error after skipping " << i - nevents
           << " file exhausted?" << endl;
      errorflag = -1;
      fileclose();
    }
    else
    {
      if (Verbosity() > 3)
      {
        cout << "Skipping evt no: " << m_Event->getEvtSequence() << endl;
      }
    }
    delete m_Event;
    m_Event = nullptr;
    nevents--;
  }
  return errorflag;
}

int Fun4AllPrdfInputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (m_SyncObject)
    {
      *mastersync = dynamic_cast<SyncObject *> (m_SyncObject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *m_SyncObject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllPrdfInputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << endl;
    cout << "Most likely your first file does not contain a SyncObject and the file" << endl;
    cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << endl;
    cout << "Change your macro and use the file opened by this input manager as first input" << endl;
    cout << "and you will be okay. Fun4All will not process the current configuration" << endl
         << endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = m_SyncObject->Different(mastersync);
  if (iret)
  {
    cout << "big problem" << endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

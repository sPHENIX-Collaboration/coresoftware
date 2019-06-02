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

#include <cstdlib>
#include <iostream>                   // for operator<<, basic_ostream, endl
#include <utility>                    // for pair

using namespace std;

Fun4AllPrdfInputManager::Fun4AllPrdfInputManager(const string &name, const string &prdfnodename, const string &topnodename)
  : Fun4AllInputManager(name, prdfnodename, topnodename)
  , segment(-999)
  , events_total(0)
  , events_thisfile(0)
  , evt(nullptr)
  , save_evt(nullptr)
  , eventiterator(nullptr)
  , syncobject(new SyncObjectv1())
  , m_PrdfNodeName(prdfnodename)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  topNode = se->topNode(TopNodeName());
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  if (!PrdfNode)
  {
    PHDataNode<Event> *newNode = new PHDataNode<Event>(evt, m_PrdfNodeName, "Event");
    topNode->addNode(newNode);
  }
  return;
}

Fun4AllPrdfInputManager::~Fun4AllPrdfInputManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete syncobject;
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
  eventiterator = new fileEventiterator(fname.c_str(), status);
  events_thisfile = 0;
  if (status)
  {
    delete eventiterator;
    eventiterator = nullptr;
    cout << PHWHERE << Name() << ": could not open file " << fname << endl;
    return -1;
  }
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(fname);
  segment = runseg.second;
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
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  if (save_evt)  // if an event was pushed back, copy saved pointer and reset save_evt pointer
  {
    evt = save_evt;
    save_evt = nullptr;
    events_thisfile--;
    events_total--;
  }
  else
  {
    evt = eventiterator->getNextEvent();
  }
  PrdfNode->setData(evt);
  if (!evt)
  {
    fileclose();
    goto readagain;
  }
  if (Verbosity() > 1)
  {
    cout << Name() << " PRDF run " << evt->getRunNumber() << ", evt no: " << evt->getEvtSequence() << endl;
  }
  events_total++;
  events_thisfile++;
  SetRunNumber(evt->getRunNumber());
  MySyncManager()->PrdfEvents(events_thisfile);
  MySyncManager()->SegmentNumber(segment);
  MySyncManager()->CurrentEvent(evt->getEvtSequence());
  syncobject->EventCounter(events_thisfile);
  syncobject->SegmentNumber(segment);
  syncobject->RunNumber(evt->getRunNumber());
  syncobject->EventNumber(evt->getEvtSequence());
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
  delete eventiterator;
  eventiterator = nullptr;
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
  PHNodeIterator iter(topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  PrdfNode->setData(nullptr);  // set pointer in Node to nullptr before deleting it
  delete evt;
  evt = nullptr;
  syncobject->Reset();
  return 0;
}

int Fun4AllPrdfInputManager::PushBackEvents(const int i)
{
  // PushBackEvents is supposedly pushing events back on the stack which works
  // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
  // A special case is when the synchronization fails and we need to only push back a single
  // event. In this case we save the evt pointer as save_evt which is used in the run method
  // instead of getting the next event.
  if (i > 0)
  {
    if (i == 1 && evt)  // check on evt pointer makes sure it is not done from the cmd line
    {
      save_evt = evt;
      return 0;
    }
    cout << PHWHERE << Name()
         << " Fun4AllPrdfInputManager cannot push back " << i << " events into file"
         << endl;
    return -1;
  }
  if (!eventiterator)
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
    evt = eventiterator->getNextEvent();
    if (!evt)
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
        cout << "Skipping evt no: " << evt->getEvtSequence() << endl;
      }
    }
    delete evt;
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
    if (syncobject) *mastersync = syncobject->clone();
  }
  else
  {
    *(*mastersync) = *syncobject;  // copy syncobject content
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
  int iret = syncobject->Different(mastersync);
  if (iret)
  {
    cout << "big problem" << endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

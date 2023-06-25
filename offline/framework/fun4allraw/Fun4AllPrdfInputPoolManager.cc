#include "Fun4AllPrdfInputPoolManager.h"

#include "SinglePrdfInput.h"

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <ffaobjects/SyncObject.h>  // for SyncObject
#include <ffaobjects/SyncObjectv1.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/phool.h>           // for PHWHERE

#include <Event/A_Event.h>
#include <Event/Event.h>
#include <Event/Eventiterator.h>  // for Eventiterator
#include <Event/fileEventiterator.h>
#include <Event/oEvent.h>

#include <cassert>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllPrdfInputPoolManager::Fun4AllPrdfInputPoolManager(const std::string &name, const std::string &prdfnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, prdfnodename, topnodename)
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
  oph = new oEvent(workmem, 4 * 1024 * 1024, 1, 1, 1);
  return;
}

Fun4AllPrdfInputPoolManager::~Fun4AllPrdfInputPoolManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete m_SyncObject;
  for (auto iter : m_PrdfInputVector)
  {
    delete iter;
  }
  for (auto pktinfoiter : m_PacketMap)
  {
    for (auto &pktiter : pktinfoiter.second.PacketVector)
    {
      delete pktiter;
    }
  }
  delete oph;
}

int Fun4AllPrdfInputPoolManager::run(const int /*nevents*/)
{
  if (m_PacketMap.size() < 5)
  {
    for (auto iter : m_PrdfInputVector)
    {
      iter->FillPool(5);
      m_RunNumber = iter->RunNumber();
    }
    SetRunNumber(m_RunNumber);
  }

  if(m_PacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
//  std::cout << "next event is " << m_PacketMap.begin()->first << std::endl;
  auto pktinfoiter = m_PacketMap.begin();
  oph->prepare_next(pktinfoiter->first, m_RunNumber);
  for (auto &pktiter : pktinfoiter->second.PacketVector)
  {
    oph->addPacket(pktiter);
  }
  m_Event = new A_Event(workmem);
  if (Verbosity() > 1)
  {
    m_Event->identify();
  }
  PHNodeIterator iter(m_topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  PrdfNode->setData(m_Event);
  for (auto &pktiter : pktinfoiter->second.PacketVector)
  {
    delete pktiter;
  }
  m_PacketMap.erase(pktinfoiter);

  return 0;
  // readagain:
  //   if (!IsOpen())
  //   {
  //     if (FileListEmpty())
  //     {
  //       if (Verbosity() > 0)
  //       {
  //         std::cout << Name() << ": No Input file open" << std::endl;
  //       }
  //       return -1;
  //     }
  //     else
  //     {
  //       if (OpenNextFile())
  //       {
  //         std::cout << Name() << ": No Input file from filelist opened" << std::endl;
  //         return -1;
  //       }
  //     }
  //   }
  //   if (Verbosity() > 3)
  //   {
  //     std::cout << "Getting Event from " << Name() << std::endl;
  //   }
  // // Fill Event combiner
  //   unsigned int watermark = m_EventCombiner.size();
  //   if (watermark < m_LowWaterMark)
  //   {
  //     for (unsigned int i = watermark; i < m_CombinerDepth; i++)
  //     {
  //       Event *evt = m_EventIterator->getNextEvent();
  //       std::cout << "Filling combiner with event " << evt->getEvtSequence() << std::endl;
  //       m_EventCombiner.insert(std::make_pair(evt->getEvtSequence(), evt));
  //     }
  //   }
  //   //  std::cout << "running event " << nevents << std::endl;
  //   PHNodeIterator iter(m_topNode);
  //   PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  //   if (m_SaveEvent)  // if an event was pushed back, copy saved pointer and reset m_SaveEvent pointer
  //   {
  //     m_Event = m_SaveEvent;
  //     m_SaveEvent = nullptr;
  //     m_EventsThisFile--;
  //     m_EventsTotal--;
  //   }
  //   else
  //   {
  //     m_Event = m_EventCombiner.begin()->second;
  //   }
  //   PrdfNode->setData(m_Event);
  //   if (!m_Event)
  //   {
  //     fileclose();
  //     goto readagain;
  //   }
  //   if (Verbosity() > 1)
  //   {
  //     std::cout << Name() << " PRDF run " << m_Event->getRunNumber() << ", evt no: " << m_Event->getEvtSequence() << std::endl;
  //   }
  //   m_EventsTotal++;
  //   m_EventsThisFile++;
  //   SetRunNumber(m_Event->getRunNumber());
  //   MySyncManager()->PrdfEvents(m_EventsThisFile);
  //   MySyncManager()->SegmentNumber(m_Segment);
  //   MySyncManager()->CurrentEvent(m_Event->getEvtSequence());
  //   m_SyncObject->EventCounter(m_EventsThisFile);
  //   m_SyncObject->SegmentNumber(m_Segment);
  //   m_SyncObject->RunNumber(m_Event->getRunNumber());
  //   m_SyncObject->EventNumber(m_Event->getEvtSequence());
  //   // check if the local SubsysReco discards this event
  //   if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
  //   {
  //     ResetEvent();
  //     goto readagain;
  //   }
  //  return 0;
}

int Fun4AllPrdfInputPoolManager::fileclose()
{
  for (auto iter : m_PrdfInputVector)
  {
    delete iter;
  }
  m_PrdfInputVector.clear();
  return 0;
}

void Fun4AllPrdfInputPoolManager::Print(const std::string &what) const
{
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllPrdfInputPoolManager::ResetEvent()
{
  PHNodeIterator iter(m_topNode);
  PHDataNode<Event> *PrdfNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_PrdfNodeName));
  PrdfNode->setData(nullptr);  // set pointer in Node to nullptr before deleting it
  delete m_Event;
  m_Event = nullptr;
  //  m_SyncObject->Reset();
  return 0;
}

int Fun4AllPrdfInputPoolManager::PushBackEvents(const int /*i*/)
{
  return 0;
  // PushBackEvents is supposedly pushing events back on the stack which works
  // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
  // A special case is when the synchronization fails and we need to only push back a single
  // event. In this case we save the m_Event pointer as m_SaveEvent which is used in the run method
  // instead of getting the next event.
  // if (i > 0)
  // {
  //   if (i == 1 && m_Event)  // check on m_Event pointer makes sure it is not done from the cmd line
  //   {
  //     m_SaveEvent = m_Event;
  //     return 0;
  //   }
  //   std::cout << PHWHERE << Name()
  //        << " Fun4AllPrdfInputPoolManager cannot push back " << i << " events into file"
  //        << std::endl;
  //   return -1;
  // }
  // if (!m_EventIterator)
  // {
  //   std::cout << PHWHERE << Name()
  //        << " no file open" << std::endl;
  //   return -1;
  // }
  // // Skipping events is implemented as
  // // pushing a negative number of events on the stack, so in order to implement
  // // the skipping of events we read -i events.
  // int nevents = -i;  // negative number of events to push back -> skip num events
  // int errorflag = 0;
  // while (nevents > 0 && !errorflag)
  // {
  //   m_Event = m_EventIterator->getNextEvent();
  //   if (!m_Event)
  //   {
  //     std::cout << "Error after skipping " << i - nevents
  //          << " file exhausted?" << std::endl;
  //     errorflag = -1;
  //     fileclose();
  //   }
  //   else
  //   {
  //     if (Verbosity() > 3)
  //     {
  //       std::cout << "Skipping evt no: " << m_Event->getEvtSequence() << std::endl;
  //     }
  //   }
  //   delete m_Event;
  //   m_Event = nullptr;
  //   nevents--;
  // }
  // return errorflag;
}

int Fun4AllPrdfInputPoolManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (m_SyncObject)
    {
      *mastersync = dynamic_cast<SyncObject *>(m_SyncObject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *m_SyncObject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllPrdfInputPoolManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    std::cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << std::endl;
    std::cout << "Most likely your first file does not contain a SyncObject and the file" << std::endl;
    std::cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << std::endl;
    std::cout << "Change your macro and use the file opened by this input manager as first input" << std::endl;
    std::cout << "and you will be okay. Fun4All will not process the current configuration" << std::endl
              << std::endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = m_SyncObject->Different(mastersync);
  if (iret)
  {
    std::cout << "big problem" << std::endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

std::string Fun4AllPrdfInputPoolManager::GetString(const std::string &what) const
{
  if (what == "PRDFNODENAME")
  {
    return m_PrdfNodeName;
  }
  return "";
}

SinglePrdfInput *Fun4AllPrdfInputPoolManager::AddPrdfInputFile(const std::string &filenam)
{
  SinglePrdfInput *prdfin = new SinglePrdfInput("PRDFIN_" + std::to_string(m_PrdfInputVector.size()), this);
  prdfin->AddFile(filenam);
  m_PrdfInputVector.push_back(prdfin);
  return m_PrdfInputVector.back();
}

SinglePrdfInput *Fun4AllPrdfInputPoolManager::AddPrdfInputList(const std::string &filenam)
{
  SinglePrdfInput *prdfin = new SinglePrdfInput("PRDFIN_" + std::to_string(m_PrdfInputVector.size()), this);
  prdfin->AddListFile(filenam);
  m_PrdfInputVector.push_back(prdfin);
  return m_PrdfInputVector.back();
}

void Fun4AllPrdfInputPoolManager::AddPacket(const int evtno, Packet *p)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding packet " << p->getIdentifier() << " to event no " << evtno << std::endl;
  }
  m_PacketMap[evtno].PacketVector.push_back(p);
}

void Fun4AllPrdfInputPoolManager::UpdateEventFoundCounter(const int evtno)
{
  m_PacketMap[evtno].EventFoundCounter++;
}

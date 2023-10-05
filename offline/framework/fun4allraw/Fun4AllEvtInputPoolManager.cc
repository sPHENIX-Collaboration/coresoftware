#include "Fun4AllEvtInputPoolManager.h"

#include "SingleInttInput.h"
#include "SingleEvtInput.h"
#include "SingleTpcInput.h"

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainerv1.h>
#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainerv1.h>
#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainerv1.h>

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
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Event/A_Event.h>
#include <Event/Event.h>
#include <Event/Eventiterator.h>  // for Eventiterator
#include <Event/fileEventiterator.h>

#include <TSystem.h>

#include <cassert>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllEvtInputPoolManager::Fun4AllEvtInputPoolManager(const std::string &name, const std::string &dstnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, dstnodename, topnodename)
  , m_SyncObject(new SyncObjectv1())
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());
  return;
}

Fun4AllEvtInputPoolManager::~Fun4AllEvtInputPoolManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete m_SyncObject;
// clear leftover event maps
  for (auto const &mapiter : m_InttRawHitMap)
  {
    for (auto intthititer :  mapiter.second.InttRawHitVector)
    {
      delete intthititer;
    }
  }
  m_InttRawHitMap.clear();

  for (auto iter : m_EvtInputVector)
  {
    delete iter;
  }
}

int Fun4AllEvtInputPoolManager::run(const int /*nevents*/)
{
  int iret = 0;
  if (m_intt_registered_flag)
  {
    iret += FillIntt();
  }
  if (m_tpc_registered_flag)
  {
    iret += FillTpc();
  }
  
  if (m_micromegas_registered_flag)
  {
    iret += FillMicromegas();
  }

  // std::cout << "size  m_InttRawHitMap: " <<  m_InttRawHitMap.size()
  // 	    << std::endl;
  return iret;
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
  //   PHDataNode<Event> *EvtNode = dynamic_cast<PHDataNode<Event> *>(iter.findFirst("PHDataNode", m_EvtNodeName));
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
  //   EvtNode->setData(m_Event);
  //   if (!m_Event)
  //   {
  //     fileclose();
  //     goto readagain;
  //   }
  //   if (Verbosity() > 1)
  //   {
  //     std::cout << Name() << " EVT run " << m_Event->getRunNumber() << ", evt no: " << m_Event->getEvtSequence() << std::endl;
  //   }
  //   m_EventsTotal++;
  //   m_EventsThisFile++;
  //   SetRunNumber(m_Event->getRunNumber());
  //   MySyncManager()->EvtEvents(m_EventsThisFile);
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

int Fun4AllEvtInputPoolManager::fileclose()
{
  for (auto iter : m_EvtInputVector)
  {
    delete iter;
  }
  m_EvtInputVector.clear();
  return 0;
}

void Fun4AllEvtInputPoolManager::Print(const std::string &what) const
{
  std::cout << "Current list of beamclks: " << std::endl;
  for (auto const &mapiter : m_PacketInfoMap)
  {
    std::cout << "clk: 0x" << std::hex << mapiter.first
              << std::dec << std::endl;
    for (auto pktiter : mapiter.second.PacketVector)
    {
      std::cout << "pktid: " << pktiter->getIdentifier() << std::endl;
    }
  }
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllEvtInputPoolManager::ResetEvent()
{
  // for (auto iter : m_EvtInputVector)
  // {
  //   iter->CleanupUsedPackets(m_CurrentBeamClock);
  // }
  //  m_SyncObject->Reset();
  return 0;
}

int Fun4AllEvtInputPoolManager::PushBackEvents(const int /*i*/)
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
  //        << " Fun4AllEvtInputPoolManager cannot push back " << i << " events into file"
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

int Fun4AllEvtInputPoolManager::GetSyncObject(SyncObject **mastersync)
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

int Fun4AllEvtInputPoolManager::SyncIt(const SyncObject *mastersync)
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

std::string Fun4AllEvtInputPoolManager::GetString(const std::string &what) const
{
  std::cout << PHWHERE << " called with " << what << " , returning empty string" << std::endl;
  return "";
}
/*
SingleEvtInput *Fun4AllEvtInputPoolManager::AddEvtInputFile(const std::string &filenam)
{
  SingleEvtInput *evtin = new SingleEvtInput("EVTIN_" + std::to_string(m_EvtInputVector.size()), this);
  evtin->AddFile(filenam);
  m_EvtInputVector.push_back(evtin);
  return m_EvtInputVector.back();
}

SingleEvtInput *Fun4AllEvtInputPoolManager::AddEvtInputList(const std::string &filenam)
{
  SingleEvtInput *evtin = new SingleEvtInput("EVTIN_" + std::to_string(m_EvtInputVector.size()), this);
  evtin->AddListFile(filenam);
  m_EvtInputVector.push_back(evtin);
  return m_EvtInputVector.back();
}
*/
void Fun4AllEvtInputPoolManager::registerStreamingInput(SingleStreamingInput *evtin, enu_subsystem system)
{
  m_EvtInputVector.push_back(evtin);
  evtin->InputManager(this);
  evtin->CreateDSTNode(m_topNode);
  switch (system)
  {
  case Fun4AllEvtInputPoolManager::MVTX:
    m_mvtx_registered_flag = true;
    break;
  case Fun4AllEvtInputPoolManager::INTT:
    m_intt_registered_flag = true;
    break;
  case Fun4AllEvtInputPoolManager::TPC:
    m_tpc_registered_flag = true;
    break;
  case Fun4AllEvtInputPoolManager::MICROMEGAS:
    m_micromegas_registered_flag = true;
    break;
  default:
    std::cout << "invalid subsystem flag " << system << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  if (Verbosity() > 3)
  {
    std::cout << "registering " << evtin->Name()
	      << " number of registered inputs: " << m_EvtInputVector.size()
	      << std::endl;
  }
}

void Fun4AllEvtInputPoolManager::AddPacket(uint64_t bclk, Packet *p)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding packet " << p->getIdentifier() << " to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_PacketInfoMap[bclk].PacketVector.push_back(p);
}

void Fun4AllEvtInputPoolManager::AddInttRawHit(uint64_t bclk, InttRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding intt hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_InttRawHitMap[bclk].InttRawHitVector.push_back(hit);
}

void Fun4AllEvtInputPoolManager::AddMicromegasRawHit(uint64_t bclk, MicromegasRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding tpc hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_MicromegasRawHitMap[bclk].MicromegasRawHitVector.push_back(hit);
}

void Fun4AllEvtInputPoolManager::AddTpcRawHit(uint64_t bclk, TpcRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding tpc hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_TpcRawHitMap[bclk].TpcRawHitVector.push_back(hit);
}

void Fun4AllEvtInputPoolManager::UpdateEventFoundCounter(const int evtno)
{
  m_PacketInfoMap[evtno].EventFoundCounter++;
}

int Fun4AllEvtInputPoolManager::FillIntt()
{
    while (m_InttRawHitMap.size() < 5) // pooling at least 5 events
    {
      unsigned int alldone = 0;
      for (auto iter : m_EvtInputVector)
      {
	alldone += iter->AllDone();
	if (Verbosity() > 0)
	{
	  std::cout << "Fun4AllEvtInputPoolManager::FillIntt - fill pool for " << iter->Name() << std::endl;
	}
	iter->FillPool();
	m_RunNumber = iter->RunNumber();
      }
      if (alldone >= m_EvtInputVector.size())
      {
	break;
      }
      SetRunNumber(m_RunNumber);
    }
    if (m_InttRawHitMap.empty())
    {
      std::cout << "we are done" << std::endl;
      return -1;
    }
    InttRawHitContainer *inttcont =  findNode::getClass<InttRawHitContainer>(m_topNode,"INTTRAWHIT");
//  std::cout << "before filling m_InttRawHitMap size: " <<  m_InttRawHitMap.size() << std::endl;
    for (auto intthititer :  m_InttRawHitMap.begin()->second.InttRawHitVector)
    {
      if (Verbosity() > 1)
      {
	intthititer->identify();
      }
      inttcont->AddHit(intthititer);
//     delete intthititer; // cleanup up done in Single Input Mgrs
    }
    for (auto iter : m_EvtInputVector)
    {
      iter->CleanupUsedPackets(m_InttRawHitMap.begin()->first);
    }
    m_InttRawHitMap.begin()->second.InttRawHitVector.clear();
    m_InttRawHitMap.erase(m_InttRawHitMap.begin());
  // std::cout << "size  m_InttRawHitMap: " <<  m_InttRawHitMap.size()
  // 	    << std::endl;
  return 0;
}

//_______________________________________________________
int Fun4AllEvtInputPoolManager::FillMicromegas()
{
  while (m_MicromegasRawHitMap.size() < 5) // pooling at least 5 events
  {
    unsigned int alldone = 0;
    for (auto iter : m_EvtInputVector)
    {
      alldone += iter->AllDone();
      if (Verbosity() > 0)
      { std::cout << "Fun4AllEvtInputPoolManager::FillMicromegas - fill pool for " << iter->Name() << std::endl; }
      iter->FillPool();
      m_RunNumber = iter->RunNumber();
    }
    
    if (alldone >= m_EvtInputVector.size())
    { break; }
    SetRunNumber(m_RunNumber);
  }
  
  if (m_MicromegasRawHitMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  
  auto container =  findNode::getClass<MicromegasRawHitContainer>(m_topNode,"MICROMEGASRAWHIT");
  for( auto hititer :  m_MicromegasRawHitMap.begin()->second.MicromegasRawHitVector)
  {
    if (Verbosity() > 1)
    {  hititer->identify(); }
    
    container->AddHit(hititer);
  }
  
  for (auto iter : m_EvtInputVector)
  { iter->CleanupUsedPackets(m_MicromegasRawHitMap.begin()->first); }
  
  m_MicromegasRawHitMap.begin()->second.MicromegasRawHitVector.clear();
  m_MicromegasRawHitMap.erase(m_MicromegasRawHitMap.begin());
  return 0;
}

int Fun4AllEvtInputPoolManager::FillTpc()
{
    while (m_TpcRawHitMap.size() < 5) // pooling at least 5 events
    {
      unsigned int alldone = 0;
      for (auto iter : m_EvtInputVector)
      {
	alldone += iter->AllDone();
	if (Verbosity() > 0)
	{
	  std::cout << "Fun4AllEvtInputPoolManager::FillTpc - fill pool for " << iter->Name() << std::endl;
	}
	iter->FillPool();
	m_RunNumber = iter->RunNumber();
      }
      if (alldone >= m_EvtInputVector.size())
      {
	break;
      }
      SetRunNumber(m_RunNumber);
    }
    if (m_TpcRawHitMap.empty())
    {
      std::cout << "we are done" << std::endl;
      return -1;
    }
    TpcRawHitContainer *tpccont =  findNode::getClass<TpcRawHitContainer>(m_topNode,"TPCRAWHIT");
//  std::cout << "before filling m_TpcRawHitMap size: " <<  m_TpcRawHitMap.size() << std::endl;
    for (auto tpchititer :  m_TpcRawHitMap.begin()->second.TpcRawHitVector)
    {
      if (Verbosity() > 1)
      {
	tpchititer->identify();
      }
      tpccont->AddHit(tpchititer);
//     delete tpchititer; // cleanup up done in Single Input Mgrs
    }
    for (auto iter : m_EvtInputVector)
    {
      iter->CleanupUsedPackets(m_TpcRawHitMap.begin()->first);
    }
    m_TpcRawHitMap.begin()->second.TpcRawHitVector.clear();
    m_TpcRawHitMap.erase(m_TpcRawHitMap.begin());
  // std::cout << "size  m_TpcRawHitMap: " <<  m_TpcRawHitMap.size()
  // 	    << std::endl;
  return 0;
}

#include "Fun4AllPrdfInputTriggerManager.h"

#include "SinglePrdfInput.h"
#include "SingleTriggerInput.h"

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <ffaobjects/SyncObject.h>    // for SyncObject
#include <ffaobjects/SyncObjectv1.h>  // for SyncObject

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cassert>
#include <climits>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllPrdfInputTriggerManager::Fun4AllPrdfInputTriggerManager(const std::string &name, const std::string &prdfnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, prdfnodename, topnodename)
  , m_SyncObject(new SyncObjectv1())
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());
}

Fun4AllPrdfInputTriggerManager::~Fun4AllPrdfInputTriggerManager()
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
  for (auto const &mapiter : m_Gl1PacketMap)
  {
    for (auto gl1packet : mapiter.second.Gl1PacketVector)
    {
      delete gl1packet;
    }
  }
  m_Gl1PacketMap.clear();
 
}

int Fun4AllPrdfInputTriggerManager::run(const int /*nevents*/)
{
  int iret = 0;
  if (m_gl1_registered_flag)  // Gl1 first to get the reference
  {
    iret += FillGl1();
  }
  if (m_mbd_registered_flag)  // Mbd next to get the reference if Gl1 is missing
  {
    iret += FillMbd();
  }
  if (m_hcal_registered_flag)  // Mbd first to get the reference
  {
    iret += FillHcal();
  }
  if (m_cemc_registered_flag)  // Mbd first to get the reference
  {
    iret += FillCemc();
  }
  if (m_zdc_registered_flag)  // Mbd first to get the reference
  {
    iret += FillZdc();
  }
  DetermineReferenceEventNumber();
  std::cout << "new ref event: " <<  m_RefEventNo << std::endl;
  MoveGl1ToNodeTree();
  MoveMbdToNodeTree();
  MoveCemcToNodeTree();
  MoveHcalToNodeTree();
  MoveZdcToNodeTree();
  MySyncManager()->CurrentEvent(m_RefEventNo);
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

int Fun4AllPrdfInputTriggerManager::fileclose()
{
  for (auto iter : m_PrdfInputVector)
  {
    delete iter;
  }
  m_PrdfInputVector.clear();
  return 0;
}

void Fun4AllPrdfInputTriggerManager::Print(const std::string &what) const
{
  //  Fun4AllInputManager::Print(what);
  if (what == "ALL" || what == "DROPPED")
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "dropped packets:" << std::endl;
    for (auto iter : m_DroppedPacketMap)
    {
      std::cout << "Packet " << iter.first << " was dropped " << iter.second << " times" << std::endl;
    }
  }
  if (what == "ALL" || what == "INPUTFILES")
  {
    std::cout << "-----------------------------" << std::endl;
    for (const auto &iter : m_Gl1InputVector)
    {
      std::cout << "Single Prdf Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
		<< " from file " << iter->FileName()
		<< std::endl;
    }
    for (const auto &iter : m_MbdInputVector)
    { 
      std::cout << "Single Prdf Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
		<< " from file " << iter->FileName()
		<< std::endl;
    }
    for (const auto &iter : m_ZdcInputVector)
    {
      std::cout << "Single Prdf Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
		<< " from file " << iter->FileName()
		<< std::endl;
    }
    for (const auto &iter : m_CemcInputVector)
    {
      std::cout << "Single Prdf Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
		<< " from file " << iter->FileName()
		<< std::endl;
    }
    for (const auto &iter : m_HcalInputVector)
    {
      std::cout << "Single Prdf Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
		<< " from file " << iter->FileName()
		<< std::endl;
    }
  }
  return;
}

int Fun4AllPrdfInputTriggerManager::ResetEvent()
{
  m_RefEventNo = std::numeric_limits<int>::min();
  return 0;
}

int Fun4AllPrdfInputTriggerManager::PushBackEvents(const int /*i*/)
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
  //        << " Fun4AllPrdfInputTriggerManager cannot push back " << i << " events into file"
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

int Fun4AllPrdfInputTriggerManager::GetSyncObject(SyncObject **mastersync)
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

int Fun4AllPrdfInputTriggerManager::SyncIt(const SyncObject *mastersync)
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

std::string Fun4AllPrdfInputTriggerManager::GetString(const std::string &what) const
{
  std::cout << PHWHERE << " called with " << what << " , returning empty string" << std::endl;
  return "";
}


void Fun4AllPrdfInputTriggerManager::registerTriggerInput(SingleTriggerInput *prdfin, InputManagerType::enu_subsystem system)
{
  prdfin->CreateDSTNode(m_topNode);
  prdfin->TriggerInputManager(this);
  switch (system)
  {
  case InputManagerType::GL1:
    m_gl1_registered_flag = true;
    m_Gl1InputVector.push_back(prdfin);
    break;
  case InputManagerType::MBD:
    m_mbd_registered_flag = true;
    m_MbdInputVector.push_back(prdfin);
    break;
  case InputManagerType::HCAL:
    m_hcal_registered_flag = true;
    m_HcalInputVector.push_back(prdfin);
    break;
  case InputManagerType::CEMC:
    m_cemc_registered_flag = true;
    m_CemcInputVector.push_back(prdfin);
    break;
  case InputManagerType::ZDC:
    m_zdc_registered_flag = true;
    m_ZdcInputVector.push_back(prdfin);
    break;
  default:
    std::cout << "invalid subsystem flag " << system << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  if (Verbosity() > 3)
  {
    std::cout << "registering " << prdfin->Name()
              << " number of registered inputs: "
              << m_Gl1InputVector.size() + m_MbdInputVector.size()
              << std::endl;
  }

  return;
}

void Fun4AllPrdfInputTriggerManager::AddBeamClock(const int evtno, const int bclk, SinglePrdfInput *prdfin)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding event " << evtno << ", clock 0x" << std::hex << bclk << std::dec
              << " snglinput: " << prdfin->Name() << std::endl;
  }
  m_ClockCounters[evtno].push_back(std::make_pair(bclk, prdfin));
}

void Fun4AllPrdfInputTriggerManager::UpdateEventFoundCounter(const int /*evtno*/)
{
//  m_PacketMap[evtno].EventFoundCounter++;
}

void Fun4AllPrdfInputTriggerManager::UpdateDroppedPacket(const int packetid)
{
  m_DroppedPacketMap[packetid]++;
}

void Fun4AllPrdfInputTriggerManager::SetReferenceClock(const int evtno, const int bclk)
{
  m_RefClockCounters[evtno] = bclk;
}

void Fun4AllPrdfInputTriggerManager::CreateBclkOffsets()
{
  if (!m_RefPrdfInput)
  {
    std::cout << PHWHERE << " No reference input manager given" << std::endl;
    exit(1);
  }
  std::map<SinglePrdfInput *, std::map<int, int>> clockcounters;
  for (const auto &iter : m_ClockCounters)
  {
    int refclock = m_RefClockCounters[iter.first];
    for (auto veciter : iter.second)
    {
      int diffclk = CalcDiffBclk(veciter.first, refclock);
      if (Verbosity() > 1)
      {
        std::cout << "diffclk for " << veciter.second->Name() << ": " << std::hex
                  << diffclk << ", clk: 0x" << veciter.first
                  << ", refclk: 0x" << refclock << std::dec << std::endl;
      }
      auto clkiter = clockcounters.find(veciter.second);
      if (clkiter == clockcounters.end())
      {
        std::map<int, int> mymap;
        clkiter = clockcounters.insert(std::make_pair(veciter.second, mymap)).first;
      }
      clkiter->second[diffclk]++;
    }
  }
  // now loop over the clock counter diffs for each input manager and find the majority vote
  for (const auto &iter : clockcounters)
  {
    int imax = -1;
    int diffmax = INT_MAX;
    for (auto initer : iter.second)
    {
      if (Verbosity() > 0)
      {
        std::cout << iter.first->Name() << " initer.second " << initer.second << std::hex
                  << " initer.first: " << initer.first << std::dec << std::endl;
      }
      if (initer.second > imax)
      {
        diffmax = initer.first;
        imax = initer.second;
      }
      m_SinglePrdfInputInfo[iter.first].bclkoffset = diffmax;
    }
  }
  for (auto iter : m_SinglePrdfInputInfo)
  {
    if (Verbosity() > 0)
    {
      std::cout << "prdf mgr " << iter.first->Name() << " clkdiff: 0x" << std::hex
                << iter.second.bclkoffset << std::dec << std::endl;
    }
  }
}

uint64_t Fun4AllPrdfInputTriggerManager::CalcDiffBclk(const uint64_t bclk1, const uint64_t bclk2)
{
  uint64_t diffclk = (bclk2 - bclk1) & 0xFFFFU;
  return diffclk;
}

void Fun4AllPrdfInputTriggerManager::DitchEvent(const int eventno)
{
  if (Verbosity() > 1)
  {
    std::cout << "Killing event " << eventno << std::endl;
  }
  return;
/*
  m_ClockCounters.erase(eventno);
  m_RefClockCounters.erase(eventno);
  auto pktinfoiter = m_PacketMap.find(eventno);
  if (pktinfoiter == m_PacketMap.end())
  {
    return;
  }
  for (auto const &pktiter : pktinfoiter->second.PacketVector)
  {
    delete pktiter;
  }
  m_PacketMap.erase(pktinfoiter);
  return;
*/
}
void Fun4AllPrdfInputTriggerManager::Resynchronize()
{
  // just load events to give us a chance to find the match
  struct LocalInfo
  {
    int clockcounter;
    int eventdiff;
  };
  for (auto iter : m_PrdfInputVector)
  {
    iter->FillPool(10);
    //    iter->FillPool(m_InitialPoolDepth);
    m_RunNumber = iter->RunNumber();
  }
  std::map<SinglePrdfInput *, LocalInfo> matchevent;
  std::vector<int> ditchevents;
  for (auto iter : m_RefClockCounters)
  {
    if (Verbosity() > 1)
    {
      std::cout << "looking for matching event " << iter.first
                << std::hex << " with clk 0x" << iter.second << std::dec << std::endl;
    }
    for (const auto &clockiter : m_ClockCounters)
    {
      if (Verbosity() > 1)
      {
        std::cout << "testing for matching with event " << clockiter.first << std::endl;
      }
      for (auto eventiter : clockiter.second)
      {
        uint64_t diffclock = CalcDiffBclk(eventiter.first, iter.second);
        if (Verbosity() > 1)
        {
          std::cout << "Event " << iter.first << " match with event " << clockiter.first
                    << " clock 0x" << std::hex << eventiter.first << ", ref clock 0x" << iter.second
                    << " diff 0x" << diffclock << std::dec
                    << " for " << eventiter.second->Name() << std::endl;
        }
        if (diffclock == m_SinglePrdfInputInfo[eventiter.second].bclkoffset)
        {
          if (Verbosity() > 1)
          {
            std::cout << "looking good for " << eventiter.second->Name() << std::endl;
          }
          matchevent[eventiter.second].clockcounter = clockiter.first;
          matchevent[eventiter.second].eventdiff = clockiter.first - iter.first;
        }
        else
        {
          if (Verbosity() > 1)
          {
            std::cout << "not so great for " << eventiter.second->Name() << std::endl;
          }
        }
      }
      if (matchevent.size() == m_SinglePrdfInputInfo.size())
      {
        if (Verbosity() > 1)
        {
          std::cout << "found all matches" << std::endl;
        }
        break;
      }
    }
    if (matchevent.size() == m_SinglePrdfInputInfo.size())
    {
      if (Verbosity() > 1)
      {
        std::cout << "found all matches" << std::endl;
      }
      break;
    }
    ditchevents.push_back(iter.first);
  }
  for (auto ievent : ditchevents)
  {
    DitchEvent(ievent);
  }
  int minoffset = std::numeric_limits<int>::max();
  for (auto matches : matchevent)
  {
    if (Verbosity() > 1)
    {
      std::cout << matches.first->Name() << " update event offset with: " << matches.second.eventdiff
                << ", current offset : " << matches.first->EventNumberOffset()
                << " would go to " << matches.first->EventNumberOffset() - matches.second.eventdiff << std::endl;
    }
    if (minoffset > matches.first->EventNumberOffset() - matches.second.eventdiff)
    {
      minoffset = matches.first->EventNumberOffset() - matches.second.eventdiff;
    }
  }
  // we cannot have negative offsets right now (this would require re-reading the previous event which is gone)
  int addoffset = 0;
  if (minoffset < 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << "minoffset < 0: " << minoffset << " this will be interesting" << std::endl;
    }
    addoffset = -minoffset;
  }
  for (auto matches : matchevent)
  {
    matches.first->EventNumberOffset(matches.first->EventNumberOffset() - matches.second.eventdiff + addoffset);
    if (Verbosity() > 1)
    {
      std::cout << matches.first->Name() << " update event offset to: " << matches.first->EventNumberOffset()
                << std::endl;
    }
  }
  ClearAllEvents();
  return;
}

void Fun4AllPrdfInputTriggerManager::ClearAllEvents()
{
  for (auto const &mapiter : m_Gl1PacketMap)
  {
    for (auto gl1packet : mapiter.second.Gl1PacketVector)
    {
      delete gl1packet;
    }
  }
  m_Gl1PacketMap.clear();
  m_ClockCounters.clear();
  m_RefClockCounters.clear();
}

int Fun4AllPrdfInputTriggerManager::FillGl1()
{
  // unsigned int alldone = 0;
  for (auto iter : m_Gl1InputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillGl1 - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool();
    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
	std::cout << PHWHERE << " Run Number mismatch, run is "
		  << m_RunNumber << ", " << iter->Name() << " reads "
		  << iter->RunNumber() << std::endl;
	std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
	Print("INPUTFILES");
	gSystem->Exit(1);
	exit(1);
      }
    }
  }
  if (m_Gl1PacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveGl1ToNodeTree()
{
    std::cout << "stashed gl1 Events: " << m_Gl1PacketMap.size() << std::endl;
  Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(m_topNode, "GL1Packet");
  //  std::cout << "before filling m_Gl1PacketMap size: " <<  m_Gl1PacketMap.size() << std::endl;
  if (! gl1packet)
  {
    return 0;
  }

  for (auto gl1hititer : m_Gl1PacketMap.begin()->second.Gl1PacketVector)
  {
    if (Verbosity() > 1)
    {
      gl1hititer->identify();
    }
    gl1packet->FillFrom(gl1hititer);
    // m_RefEventNo = gl1hititer->getEvtSequence();
    // gl1packet->setEvtSequence(m_RefEventNo);
  }
  for (auto iter : m_Gl1InputVector)
  {
    iter->CleanupUsedPackets(m_Gl1PacketMap.begin()->first);
  }
  m_Gl1PacketMap.begin()->second.Gl1PacketVector.clear();
  m_Gl1PacketMap.erase(m_Gl1PacketMap.begin());
  // std::cout << "size  m_Gl1PacketMap: " <<  m_Gl1PacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddGl1Packet(int eventno, Gl1Packet *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding gl1 hit to eventno: "
              << eventno << std::endl;
  }
  m_Gl1PacketMap[eventno].Gl1PacketVector.push_back(pkt);
  return;
}

int Fun4AllPrdfInputTriggerManager::FillMbd()
{
  // unsigned int alldone = 0;
  for (auto iter : m_MbdInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillMbd - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool();
    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
	std::cout << PHWHERE << " Run Number mismatch, run is "
		 << m_RunNumber << ", " << iter->Name() << " reads "
		  << iter->RunNumber() << std::endl;
	std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
	Print("INPUTFILES");
	gSystem->Exit(1);
	exit(1);
      }
    }
  }
  if (m_MbdPacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveMbdToNodeTree()
{
  //std::cout << "stashed mbd Events: " << m_MbdPacketMap.size() << std::endl;
  CaloPacketContainer *mbd = findNode::getClass<CaloPacketContainer>(m_topNode, "MBDPackets");
  if (! mbd)
  {
    return 0;
  }
  //  std::cout << "before filling m_MbdPacketMap size: " <<  m_MbdPacketMap.size() << std::endl;
  for (auto mbdhititer : m_MbdPacketMap.begin()->second.MbdPacketVector)
  {
    if (Verbosity() > 1)
    {
      mbdhititer->identify();
    }
    mbd->AddPacket(mbdhititer);
  }
  for (auto iter : m_MbdInputVector)
  {
    iter->CleanupUsedPackets(m_MbdPacketMap.begin()->first);
  }
  m_MbdPacketMap.begin()->second.MbdPacketVector.clear();
  m_MbdPacketMap.erase(m_MbdPacketMap.begin());
  // std::cout << "size  m_MbdPacketMap: " <<  m_MbdPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddMbdPacket(int eventno, CaloPacket *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding mbd hit to eventno: "
              << eventno << std::endl;
  }
  m_MbdPacketMap[eventno].MbdPacketVector.push_back(pkt);
  return;
}

int Fun4AllPrdfInputTriggerManager::FillHcal()
{
  // unsigned int alldone = 0;
  for (auto iter : m_HcalInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillHcal - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool();
    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
	std::cout << PHWHERE << " Run Number mismatch, run is "
		 << m_RunNumber << ", " << iter->Name() << " reads "
		  << iter->RunNumber() << std::endl;
	std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
	Print("INPUTFILES");
	gSystem->Exit(1);
	exit(1);
      }
    }
  }
  if (m_HcalPacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveHcalToNodeTree()
{
  std::cout << "stashed hcal Events: " << m_HcalPacketMap.size() << std::endl;
  CaloPacketContainer *hcal = findNode::getClass<CaloPacketContainer>(m_topNode, "HCALPackets");
  if (!hcal)
  {
    return 0;
  }
  //  std::cout << "before filling m_HcalPacketMap size: " <<  m_HcalPacketMap.size() << std::endl;
  for (auto hcalhititer : m_HcalPacketMap.begin()->second.HcalPacketVector)
  {
    if (Verbosity() > 1)
    {
      hcalhititer->identify();
    }
    hcal->AddPacket(hcalhititer);
  }
  for (auto iter : m_HcalInputVector)
  {
    iter->CleanupUsedPackets(m_HcalPacketMap.begin()->first);
  }
  m_HcalPacketMap.begin()->second.HcalPacketVector.clear();
  m_HcalPacketMap.erase(m_HcalPacketMap.begin());
  // std::cout << "size  m_HcalPacketMap: " <<  m_HcalPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddHcalPacket(int eventno, CaloPacket *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding hcal packet " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  m_HcalPacketMap[eventno].HcalPacketVector.push_back(pkt);
  return;
}

int Fun4AllPrdfInputTriggerManager::FillCemc()
{
  // unsigned int alldone = 0;
  for (auto iter : m_CemcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillCemc - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool();
    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
	std::cout << PHWHERE << " Run Number mismatch, run is "
		 << m_RunNumber << ", " << iter->Name() << " reads "
		  << iter->RunNumber() << std::endl;
	std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
	Print("INPUTFILES");
	gSystem->Exit(1);
	exit(1);
      }
    }
  }
  if (m_CemcPacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveCemcToNodeTree()
{
  std::cout << "stashed cemc Events: " << m_CemcPacketMap.size() << std::endl;
  CaloPacketContainer *cemc = findNode::getClass<CaloPacketContainer>(m_topNode, "CEMCPackets");
  if (! cemc)
  {
    return 0;
  }
  //  std::cout << "before filling m_CemcPacketMap size: " <<  m_CemcPacketMap.size() << std::endl;
  for (auto cemchititer : m_CemcPacketMap.begin()->second.CemcPacketVector)
  {
    if (Verbosity() > 1)
    {
      cemchititer->identify();
    }
    cemc->AddPacket(cemchititer);
  }
  for (auto iter : m_CemcInputVector)
  {
    iter->CleanupUsedPackets(m_CemcPacketMap.begin()->first);
  }
  m_CemcPacketMap.begin()->second.CemcPacketVector.clear();
  m_CemcPacketMap.erase(m_CemcPacketMap.begin());
  // std::cout << "size  m_CemcPacketMap: " <<  m_CemcPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddCemcPacket(int eventno, CaloPacket *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding cemc packet " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  m_CemcPacketMap[eventno].CemcPacketVector.push_back(pkt);
  return;
}

int Fun4AllPrdfInputTriggerManager::FillZdc()
{
  // unsigned int alldone = 0;
  for (auto iter : m_ZdcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillZdc - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool();
    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
	std::cout << PHWHERE << " Run Number mismatch, run is "
		 << m_RunNumber << ", " << iter->Name() << " reads "
		  << iter->RunNumber() << std::endl;
	std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
	Print("INPUTFILES");
	gSystem->Exit(1);
	exit(1);
      }
    }
  }
  if (m_ZdcPacketMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveZdcToNodeTree()
{
  std::cout << "stashed zdc Events: " << m_ZdcPacketMap.size() << std::endl;
  CaloPacketContainer *zdc = findNode::getClass<CaloPacketContainer>(m_topNode, "ZDCPackets");
  if (!zdc)
  {
    return 0;
  }
  //  std::cout << "before filling m_ZdcPacketMap size: " <<  m_ZdcPacketMap.size() << std::endl;
  for (auto zdchititer : m_ZdcPacketMap.begin()->second.ZdcPacketVector)
  {
    if (Verbosity() > 1)
    {
      zdchititer->identify();
    }
    zdc->AddPacket(zdchititer);
  }
  for (auto iter : m_ZdcInputVector)
  {
    iter->CleanupUsedPackets(m_ZdcPacketMap.begin()->first);
  }
  m_ZdcPacketMap.begin()->second.ZdcPacketVector.clear();
  m_ZdcPacketMap.erase(m_ZdcPacketMap.begin());
  // std::cout << "size  m_ZdcPacketMap: " <<  m_ZdcPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddZdcPacket(int eventno, CaloPacket *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding zdc packet " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  m_ZdcPacketMap[eventno].ZdcPacketVector.push_back(pkt);
  return;
}

void Fun4AllPrdfInputTriggerManager::DetermineReferenceEventNumber()
{
  if (!m_Gl1PacketMap.empty())
  {
    m_RefEventNo = m_Gl1PacketMap.begin()->first;
  }
  else if (!m_MbdPacketMap.empty())
  {
    m_RefEventNo = m_MbdPacketMap.begin()->first;
  }
  else if (!m_HcalPacketMap.empty())
  {
    m_RefEventNo = m_HcalPacketMap.begin()->first;
  }
  else if (!m_CemcPacketMap.empty())
  {
    m_RefEventNo = m_CemcPacketMap.begin()->first;
  }
  else if (!m_ZdcPacketMap.empty())
  {
    m_RefEventNo = m_ZdcPacketMap.begin()->first;
  }
    return;
}

#include "Fun4AllPrdfInputTriggerManager.h"

#include "SinglePrdfInput.h"
#include "SingleTriggerInput.h"

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <ffaobjects/SyncObject.h>    // for SyncObject
#include <ffaobjects/SyncObjectv1.h>  // for SyncObject

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/LL1Packet.h>
#include <ffarawobjects/LL1PacketContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <algorithm>  // for std::search
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
  for (auto iter : m_TriggerInputVector)
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << " deleting " << iter->Name() << std::endl;
    }
    delete iter;
  }
}

int Fun4AllPrdfInputTriggerManager::run(const int /*nevents*/)
{
tryagain:
  int iret = 0;
  if (m_gl1_registered_flag)  // Gl1 first to get the reference
  {
    iret += FillGl1(m_PoolDepth);
  }
  if (m_mbd_registered_flag)  // Mbd next to get the reference if Gl1 is missing
  {
    iret += FillMbd(m_PoolDepth);
  }
  if (m_hcal_registered_flag)  // Mbd first to get the reference
  {
    iret += FillHcal(m_PoolDepth);
  }
  if (m_cemc_registered_flag)  // Mbd first to get the reference
  {
    iret += FillCemc(m_PoolDepth);
  }
  if (m_zdc_registered_flag)  // Mbd first to get the reference
  {
    iret += FillZdc(m_PoolDepth);
  }
  if (m_ll1_registered_flag)  // LL1 next to get the reference if Gl1 is missing
  {
    iret += FillLL1(m_PoolDepth);
  }
  if (iret)
  {
    return -1;
  }
  m_PoolDepth = m_DefaultPoolDepth;
  DetermineReferenceEventNumber();
  if (Verbosity() > 0)
  {
    std::cout << "new ref event: " << m_RefEventNo << std::endl;
  }
  if (m_resync_flag)
  {
    //    Print("CEMCMAP");
    ClockDiffFill();
    if (ClockDiffCheck())
    {
      // this is not used yet - ClockDiffCheck() always returns zero
      // NOLINTNEXTLINE(hicpp-avoid-goto)
      goto tryagain;
    }
    //    Print("CEMCMAP");
  }
  MoveGl1ToNodeTree();
  MoveMbdToNodeTree();
  MoveCemcToNodeTree();
  MoveHcalToNodeTree();
  MoveLL1ToNodeTree();
  // do not switch the order of zdc and sepd, they use a common input manager
  // and the cleanup is done in MoveSEpdToNodeTree, if the MoveZdcToNodeTree is
  // called after that it will segfault
  MoveZdcToNodeTree();
  MoveSEpdToNodeTree();
  MySyncManager()->CurrentEvent(m_RefEventNo);
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

int Fun4AllPrdfInputTriggerManager::fileclose()
{
  for (auto iter : m_TriggerInputVector)
  {
    delete iter;
  }
  m_TriggerInputVector.clear();
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
  if (what == "CEMCMAP")
  {
    std::cout << "Printing CEMCMAP" << std::endl;
    for (auto iter : m_CemcPacketMap)
    {
      std::cout << "event " << iter.first << std::endl;
      for (auto packet_iter : iter.second.CaloSinglePacketMap)
      {
        int packet_id = packet_iter.first;
        auto bcoiter = iter.second.BcoDiffMap.find(packet_id);
        uint64_t bcodiff = 0x0;
        if (bcoiter != iter.second.BcoDiffMap.end())
        {
          bcodiff = bcoiter->second;
          std::cout << "Packet " << packet_id << ", bco: 0x" << std::hex << packet_iter.second->getBCO()
                    << " bcodiff: 0x" << bcodiff << std::dec << std::endl;
        }
        else
        {
          std::cout << "Packet " << packet_id << ", bco: 0x" << std::hex << packet_iter.second->getBCO()
                    << " has no bcodiff" << std::dec << std::endl;
        }
      }
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
  case InputManagerType::LL1:
    m_ll1_registered_flag = true;
    m_LL1InputVector.push_back(prdfin);
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
  m_TriggerInputVector.push_back(prdfin);
  // this is for convenience - we need to loop over all input managers except for the GL1
  if (system != InputManagerType::GL1)
  {
    m_NoGl1InputVector.push_back(prdfin);
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

void Fun4AllPrdfInputTriggerManager::ClearAllEvents(const int eventno)
{
  for (auto &mapiter : m_Gl1PacketMap)
  {
    // for (auto &gl1packet : mapiter.second.Gl1SinglePacketMap)
    // {
    //   delete gl1packet.second;
    // }
    mapiter.second.Gl1SinglePacketMap.clear();
    mapiter.second.BcoDiffMap.clear();
  }
  m_Gl1PacketMap.clear();

  for (auto &mapiter : m_MbdPacketMap)
  {
    // for (auto &mbdpacket : mapiter.second.CaloSinglePacketMap)
    // {
    //   delete mbdpacket.second;
    // }
    mapiter.second.CaloSinglePacketMap.clear();
    mapiter.second.BcoDiffMap.clear();
  }
  m_MbdPacketMap.clear();

  for (auto &mapiter : m_SEpdPacketMap)
  {
    // for (auto &sepdpacket : mapiter.second.CaloSinglePacketMap)
    // {
    //   delete sepdpacket.second;
    // }
    mapiter.second.CaloSinglePacketMap.clear();
    mapiter.second.BcoDiffMap.clear();
  }
  m_SEpdPacketMap.clear();

  for (auto &mapiter : m_ZdcPacketMap)
  {
    // for (auto &zdcpacket : mapiter.second.CaloSinglePacketMap)
    // {
    //   delete zdcpacket.second;
    // }
    mapiter.second.CaloSinglePacketMap.clear();
    mapiter.second.BcoDiffMap.clear();
  }
  m_ZdcPacketMap.clear();

  m_ClockCounters.clear();
  m_RefClockCounters.clear();

  for (auto iter : m_Gl1InputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
  for (auto iter : m_MbdInputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
  for (auto iter : m_HcalInputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
  for (auto iter : m_CemcInputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
  for (auto iter : m_LL1InputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
  for (auto iter : m_ZdcInputVector)
  {
    iter->CleanupUsedPackets(eventno);
  }
}

int Fun4AllPrdfInputTriggerManager::FillGl1(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_Gl1InputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillGl1 - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
    std::cout << "GL1 event stack is empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveGl1ToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed gl1 Events: " << m_Gl1PacketMap.size() << std::endl;
  }
  Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(m_topNode, "GL1Packet");
  //  std::cout << "before filling m_Gl1PacketMap size: " <<  m_Gl1PacketMap.size() << std::endl;
  if (!gl1packet)
  {
    return 0;
  }

  for (auto &gl1hititer : m_Gl1PacketMap.begin()->second.Gl1SinglePacketMap)
  {
    if (Verbosity() > 1)
    {
      gl1hititer.second->identify();
    }
    gl1packet->FillFrom(gl1hititer.second);
    // m_RefEventNo = gl1hititer->getEvtSequence();
    // gl1packet->setEvtSequence(m_RefEventNo);
  }
  for (auto iter : m_Gl1InputVector)
  {
    if (Verbosity() > 1)
    {
      std::cout << "GL1: cleaning out unused packets for event " << m_Gl1PacketMap.begin()->first << std::endl;
    }
    iter->CleanupUsedPackets(m_Gl1PacketMap.begin()->first);
  }
  m_Gl1PacketMap.begin()->second.Gl1SinglePacketMap.clear();
  if (Verbosity() > 1)
  {
    std::cout << "clearing bco diff map from " << std::endl;
    for (auto iter : m_Gl1PacketMap.begin()->second.BcoDiffMap)
    {
      std::cout << "Packet " << iter.first << " bco: 0x" << std::hex << iter.second << std::dec << std::endl;
    }
  }
  m_Gl1PacketMap.begin()->second.BcoDiffMap.clear();
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
  auto &iter = m_Gl1PacketMap[eventno];
  iter.Gl1SinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  return;
}

int Fun4AllPrdfInputTriggerManager::FillMbd(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_MbdInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillMbd - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
    std::cout << "MBD event stack is empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveMbdToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed mbd Events: " << m_MbdPacketMap.size() << std::endl;
  }
  CaloPacketContainer *mbd = findNode::getClass<CaloPacketContainer>(m_topNode, "MBDPackets");
  if (!mbd)
  {
    return 0;
  }
  //  std::cout << "before filling m_MbdPacketMap size: " <<  m_MbdPacketMap.size() << std::endl;
  mbd->setEvtSequence(m_RefEventNo);
  for (auto mbdhititer : m_MbdPacketMap.begin()->second.CaloSinglePacketMap)
  {
    if (m_MbdPacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 1)
      {
        mbdhititer.second->identify();
      }
      mbd->AddPacket(mbdhititer.second);
    }
  }
  for (auto iter : m_MbdInputVector)
  {
    iter->CleanupUsedPackets(m_MbdPacketMap.begin()->first);
  }
  m_MbdPacketMap.begin()->second.CaloSinglePacketMap.clear();
  m_MbdPacketMap.begin()->second.BcoDiffMap.clear();
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
  m_MbdPacketMap[eventno].CaloSinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  return;
}

int Fun4AllPrdfInputTriggerManager::FillHcal(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_HcalInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillHcal - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
    std::cout << "Hcal event stack is empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveHcalToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed hcal Events: " << m_HcalPacketMap.size() << std::endl;
  }
  CaloPacketContainer *hcal = findNode::getClass<CaloPacketContainer>(m_topNode, "HCALPackets");
  if (!hcal)
  {
    return 0;
  }
  //  std::cout << "before filling m_HcalPacketMap size: " <<  m_HcalPacketMap.size() << std::endl;
  hcal->setEvtSequence(m_RefEventNo);
  for (auto hcalhititer : m_HcalPacketMap.begin()->second.CaloSinglePacketMap)
  {
    if (m_HcalPacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 1)
      {
        hcalhititer.second->identify();
      }
      hcal->AddPacket(hcalhititer.second);
    }
  }
  for (auto iter : m_HcalInputVector)
  {
    iter->CleanupUsedPackets(m_HcalPacketMap.begin()->first);
  }
  m_HcalPacketMap.begin()->second.CaloSinglePacketMap.clear();
  m_HcalPacketMap.begin()->second.BcoDiffMap.clear();
  m_HcalPacketMap.erase(m_HcalPacketMap.begin());
  // std::cout << "size  m_HcalPacketMap: " <<  m_HcalPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddHcalPacket(int eventno, CaloPacket *pkt)
{
  if (pkt == nullptr)
  {
    std::cout << PHWHERE << " got null ptr to add packet, not doing this" << std::endl;
    return;
  }
  if (Verbosity() > 1)
  {
    std::cout << "Adding hcal packet " << pkt->getIdentifier() << " from event " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  auto ret = m_HcalPacketMap[eventno].CaloSinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  if (ret.second)
  {
    if (Verbosity() > 1)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " was successful" << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 3)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " failed - duplicate?" << std::endl;
    }
  }
  return;
}

int Fun4AllPrdfInputTriggerManager::FillCemc(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_CemcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillCemc - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
    std::cout << "Cemc event stack is empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveCemcToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed cemc Events: " << m_CemcPacketMap.size() << std::endl;
  }
  CaloPacketContainer *cemc = findNode::getClass<CaloPacketContainer>(m_topNode, "CEMCPackets");
  if (!cemc)
  {
    return 0;
  }
  //  std::cout << "before filling m_CemcPacketMap size: " <<  m_CemcPacketMap.size() << std::endl;
  cemc->setEvtSequence(m_RefEventNo);
  if (Verbosity() > 1)
  {
    if (m_CemcPacketMap.begin()->second.CaloSinglePacketMap.empty())
    {
      std::cout << "Event " << m_RefEventNo << " is missing from CEMC" << std::endl;
    }
  }
  for (auto cemchititer : m_CemcPacketMap.begin()->second.CaloSinglePacketMap)
  {
    if (m_CemcPacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 21)
      {
        cemchititer.second->identify();
      }
      cemc->AddPacket(cemchititer.second);
    }
  }
  for (auto iter : m_CemcInputVector)
  {
    iter->CleanupUsedPackets(m_CemcPacketMap.begin()->first);
  }
  m_CemcPacketMap.begin()->second.CaloSinglePacketMap.clear();
  m_CemcPacketMap.begin()->second.BcoDiffMap.clear();
  m_CemcPacketMap.erase(m_CemcPacketMap.begin());
  // std::cout << "size  m_CemcPacketMap: " <<  m_CemcPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddCemcPacket(int eventno, CaloPacket *pkt)
{
  if (pkt == nullptr)
  {
    std::cout << PHWHERE << " got null ptr to add packet, not doing this" << std::endl;
    return;
  }
  if (Verbosity() > 1)
  {
    std::cout << "Adding cemc packet " << pkt->getIdentifier() << " from event " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  auto ret = m_CemcPacketMap[eventno].CaloSinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  if (ret.second)
  {
    if (Verbosity() > 1)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " was successful" << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 3)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " failed - duplicate?" << std::endl;
    }
  }

  //  std::cout << "Cemc packet map size: " << m_CemcPacketMap.size() << std::endl;
  return;
}

int Fun4AllPrdfInputTriggerManager::FillLL1(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_LL1InputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillLL1 - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
  if (m_LL1PacketMap.empty())
  {
    std::cout << "LL1 event stack empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveLL1ToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed ll1 Events: " << m_LL1PacketMap.size() << std::endl;
  }
  LL1PacketContainer *ll1 = findNode::getClass<LL1PacketContainer>(m_topNode, "LL1Packets");
  if (!ll1)
  {
    return 0;
  }
  //  std::cout << "before filling m_LL1PacketMap size: " <<  m_LL1PacketMap.size() << std::endl;
  ll1->setEvtSequence(m_RefEventNo);
  for (auto ll1hititer : m_LL1PacketMap.begin()->second.LL1SinglePacketMap)
  {
    if (m_LL1PacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 1)
      {
        ll1hititer.second->identify();
      }
      ll1->AddPacket(ll1hititer.second);
    }
  }
  for (auto iter : m_LL1InputVector)
  {
    iter->CleanupUsedPackets(m_LL1PacketMap.begin()->first);
  }
  m_LL1PacketMap.begin()->second.LL1SinglePacketMap.clear();
  m_LL1PacketMap.begin()->second.BcoDiffMap.clear();
  m_LL1PacketMap.erase(m_LL1PacketMap.begin());
  // std::cout << "size  m_LL1PacketMap: " <<  m_LL1PacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddLL1Packet(int eventno, LL1Packet *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "AddLL1Packet: Adding ll1 packet " << pkt->getIdentifier()
              << " for event " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  m_LL1PacketMap[eventno].LL1SinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  return;
}

int Fun4AllPrdfInputTriggerManager::FillZdc(const unsigned int nEvents)
{
  // unsigned int alldone = 0;
  for (auto iter : m_ZdcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllTriggerInputManager::FillZdc - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(nEvents);
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
    std::cout << "Zdc/Sepd event stack is empty, we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::MoveZdcToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed zdc Events: " << m_ZdcPacketMap.size() << std::endl;
  }
  CaloPacketContainer *zdc = findNode::getClass<CaloPacketContainer>(m_topNode, "ZDCPackets");
  if (!zdc)
  {
    return 0;
  }
  //  std::cout << "before filling m_ZdcPacketMap size: " <<  m_ZdcPacketMap.size() << std::endl;
  zdc->setEvtSequence(m_RefEventNo);
  for (auto zdchititer : m_ZdcPacketMap.begin()->second.CaloSinglePacketMap)
  {
    if (m_ZdcPacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 2)
      {
        std::cout << "event at m_ZdcPacketMap.begin(): " << m_ZdcPacketMap.begin()->first << std::endl;
        if (Verbosity() > 10)
        {
          zdchititer.second->identify();
        }
      }
      zdc->AddPacket(zdchititer.second);
    }
  }
  // Since the ZDC and sEPD are in the same file using the same input manager
  // clean up zdc and sepd together in MoveSEpdToNodeTree()

  /*
    for (auto iter : m_ZdcInputVector)
    {
      iter->CleanupUsedPackets(m_ZdcPacketMap.begin()->first);
    }
    m_ZdcPacketMap.begin()->second.ZdcPacketVector.clear();
    m_ZdcPacketMap.erase(m_ZdcPacketMap.begin());
  */
  // std::cout << "size  m_ZdcPacketMap: " <<  m_ZdcPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddZdcPacket(int eventno, CaloPacket *pkt)
{
  if (pkt == nullptr)
  {
    std::cout << PHWHERE << " got null ptr to add packet, not doing this" << std::endl;
    return;
  }
  if (Verbosity() > 1)
  {
    std::cout << "AddZdcPacket: Adding zdc packet " << pkt->getIdentifier()
              << " from event " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  auto ret = m_ZdcPacketMap[eventno].CaloSinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  if (ret.second)
  {
    if (Verbosity() > 1)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " was successful" << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 3)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " failed - duplicate?" << std::endl;
    }
  }
  return;
}

int Fun4AllPrdfInputTriggerManager::MoveSEpdToNodeTree()
{
  if (Verbosity() > 1)
  {
    std::cout << "stashed sepd Events: " << m_SEpdPacketMap.size() << std::endl;
  }
  CaloPacketContainer *sepd = findNode::getClass<CaloPacketContainer>(m_topNode, "SEPDPackets");
  if (!sepd)
  {
    return 0;
  }
  // std::cout << "before filling m_SEpdPacketMap size: " <<  m_SEpdPacketMap.size() << std::endl;
  sepd->setEvtSequence(m_RefEventNo);
  for (auto sepdhititer : m_SEpdPacketMap.begin()->second.CaloSinglePacketMap)
  {
    if (m_SEpdPacketMap.begin()->first == m_RefEventNo)
    {
      if (Verbosity() > 2)
      {
        std::cout << "event at m_SEpdPacketMap.begin(): " << m_SEpdPacketMap.begin()->first << std::endl;
        if (Verbosity() > 10)
        {
          sepdhititer.second->identify();
        }
      }
      sepd->AddPacket(sepdhititer.second);
    }
  }
  // Since the ZDC and sEPD are in the same file using the same input manager
  // clean up zdc and sepd here
  for (auto iter : m_ZdcInputVector)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Cleaning event no from zdc inputmgr " << m_RefEventNo << std::endl;
    }
    iter->CleanupUsedPackets(m_RefEventNo);
  }
  if (m_ZdcPacketMap.begin()->first <= m_RefEventNo)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Erasing event no " << m_ZdcPacketMap.begin()->first << " from zdc pktmap" << std::endl;
    }
    m_ZdcPacketMap.begin()->second.CaloSinglePacketMap.clear();
    m_ZdcPacketMap.begin()->second.BcoDiffMap.clear();
    m_ZdcPacketMap.erase(m_ZdcPacketMap.begin());
  }
  if (m_SEpdPacketMap.begin()->first <= m_RefEventNo)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Erasing event no " << m_SEpdPacketMap.begin()->first << " from sepd pktmap" << std::endl;
    }
    m_SEpdPacketMap.begin()->second.CaloSinglePacketMap.clear();
    m_SEpdPacketMap.begin()->second.BcoDiffMap.clear();
    m_SEpdPacketMap.erase(m_SEpdPacketMap.begin());
  }
  // std::cout << "size  m_SEpdPacketMap: " <<  m_SEpdPacketMap.size()
  // 	    << std::endl;
  return 0;
}

void Fun4AllPrdfInputTriggerManager::AddSEpdPacket(int eventno, CaloPacket *pkt)
{
  if (Verbosity() > 1)
  {
    std::cout << "AddSEpdPacket: Adding sepd packet " << pkt->getIdentifier()
              << " for event " << pkt->getEvtSequence() << " to eventno: "
              << eventno << std::endl;
  }
  //  auto &iter = m_SEpdPacketMap[eventno];
  auto ret = m_SEpdPacketMap[eventno].CaloSinglePacketMap.insert(std::make_pair(pkt->getIdentifier(), pkt));
  if (ret.second)
  {
    if (Verbosity() > 1)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " was successful" << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 3)
    {
      std::cout << "inserting packet " << pkt->getIdentifier() << " for event " << pkt->getEvtSequence()
                << " failed - duplicate?" << std::endl;
    }
  }
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
  else if (!m_LL1PacketMap.empty())
  {
    m_RefEventNo = m_LL1PacketMap.begin()->first;
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

void Fun4AllPrdfInputTriggerManager::ClockDiffFill()
{
  // std::vector<uint64_t> BcoDiffRef;
  if (!m_Gl1PacketMap.empty())
  {
    // This is a more sophisticated way to clean the Haystack but it is sensitive to repeated bco's
    // Lets go with the simpler one by just refilling the HayStack with the stored bco diffs
    // auto gl1hititer = m_Gl1PacketMap.begin();
    // if (! gl1hititer->second.BcoDiffMap.empty()) // the very first event does not have a bco diff
    // {
    // 	auto bcoiter = gl1hititer->second.BcoDiffMap.begin();
    // 	std::cout << PHWHERE << " Haystack cleanup, bco diff for packet "
    // 		  << bcoiter->first << " is 0x" << std::hex << bcoiter->second << std::dec << std::endl;

    // 	auto haystackiter = std::find(m_HayStack.begin(),m_HayStack.end(),bcoiter->second);
    // 	if (haystackiter !=  m_HayStack.end())
    // 	{
    // 	  std::cout << PHWHERE << "found " <<  std::hex << bcoiter->second << " in haystack" << std::dec << std::endl;
    // 	}
    // 	while(*m_HayStack.begin() != bcoiter->second)
    // 	{

    // 	  std::cout << "erasing " << std::hex << *m_HayStack.begin() << std::dec << std::endl;
    // 	  m_HayStack.erase(m_HayStack.begin());
    // 	}
    // 	std::cout << "begin of haystack " << std::hex << *m_HayStack.begin() << std::dec << std::endl;
    // 	std::cout << PHWHERE << "after cleaning haystack size " << m_HayStack.size() << std::endl;
    // 	for (auto iter : m_HayStack)
    // 	{
    // 	  std::cout << "after cleaning haystack: 0x" << std::hex << iter << std::dec << std::endl;
    // 	}
    // }

    // Just clear and refill the haystack with the stored bc diffs
    // we only have one GL1 packet, so we can just take the first packet
    m_HayStack.clear();
    for (auto &gl1hititer : m_Gl1PacketMap)
    {
      if (!gl1hititer.second.BcoDiffMap.empty())  // the very first event does not have a bco diff
      {
        m_HayStack.push_back(gl1hititer.second.BcoDiffMap.begin()->second);
      }
    }
    for (auto gl1hititer = m_Gl1PacketMap.begin(); gl1hititer != m_Gl1PacketMap.end(); ++gl1hititer)
    {
      //      std::cout << "current gl1 event: " <<  gl1hititer->first << std::endl;
      // this is for the very first event, only then BcoDiffMap is empty
      // we need an entry in the haystack for every event - otherwise the counting gets really hard
      if (gl1hititer->second.BcoDiffMap.empty())
      {
        for (auto &pktiter : gl1hititer->second.Gl1SinglePacketMap)
        {
          m_HayStack.push_back(0x0);
          gl1hititer->second.BcoDiffMap[pktiter.first] = 0x0;  // this is likely not needed
        }
      }
      auto nextIt = std::next(gl1hititer);
      if (nextIt != m_Gl1PacketMap.end())
      {
        //	std::cout << "size of bcomap: " << nextIt->second.BcoDiffMap.size() << std::endl;
        if (!nextIt->second.BcoDiffMap.empty())
        {
          continue;
        }
        //	std::cout << "next gl1 event: " <<  nextIt->first << std::endl;
        for (auto &pktiter : gl1hititer->second.Gl1SinglePacketMap)
        {
          uint64_t prev_bco = pktiter.second->getBCO();
          int prev_packetid = pktiter.first;
          auto currpkt = nextIt->second.Gl1SinglePacketMap.find(prev_packetid);
          if (currpkt != nextIt->second.Gl1SinglePacketMap.end())
          {
            uint64_t curr_bco = currpkt->second->getBCO();
            uint64_t diffbco = curr_bco - prev_bco;
            nextIt->second.BcoDiffMap[prev_packetid] = diffbco;
            if (Verbosity() > 11)
            {
              std::cout << "packet " << prev_packetid << ", prev_bco 0x: " << std::hex
                        << prev_bco << ", curr_bco: 0x" << curr_bco << ", diff: 0x"
                        << diffbco << std::dec << std::endl;
              std::cout << "Pushing 0x" << std::hex << diffbco << " into haystack" << std::dec << std::endl;
            }
            m_HayStack.push_back(diffbco);
          }
        }
      }
    }
  }
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "haystack size " << m_HayStack.size() << std::endl;
    for (auto iter : m_HayStack)
    {
      std::cout << "haystack: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
  if (!m_MbdPacketMap.empty())
  {
    FillNeedle(m_MbdPacketMap.begin(), m_MbdPacketMap.end(), "mbd");
  }
  if (!m_LL1PacketMap.empty())
  {
    FillNeedleLL1(m_LL1PacketMap.begin(), m_LL1PacketMap.end(), "ll1");
  }
  if (!m_HcalPacketMap.empty())
  {
    FillNeedle(m_HcalPacketMap.begin(), m_HcalPacketMap.end(), "hcal");
  }
  if (!m_CemcPacketMap.empty())
  {
    FillNeedle(m_CemcPacketMap.begin(), m_CemcPacketMap.end(), "cemc");
  }
  if (!m_ZdcPacketMap.empty())
  {
    FillNeedle(m_ZdcPacketMap.begin(), m_ZdcPacketMap.end(), "zdc");
  }
  if (!m_SEpdPacketMap.empty())
  {
    FillNeedle(m_SEpdPacketMap.begin(), m_SEpdPacketMap.end(), "sepd");
  }
  return;
}

int Fun4AllPrdfInputTriggerManager::ClockDiffCheck()
{
  std::map<int, int> eventoffset;
  static unsigned int count = 0;
  count++;
  for (auto &iter : m_NeedleMap)
  {
    std::vector needle = iter.second;
    if (count < 2 && m_FEMClockPackets.find(iter.first) != m_FEMClockPackets.end())
    {
      if (Verbosity() > 1)
      {
        std::cout << "Packet with FEM clock issue, not doing needle matching for packet " << iter.first << std::endl;
      }
      continue;
    }
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "Initial HayStack/Needle: " << iter.first
                << " HayStack size: " << m_HayStack.size() << " Needle size: " << needle.size() << std::endl;
      for (auto &hayiter : m_HayStack)
      {
        std::cout << "haystack: 0x" << std::hex << hayiter << std::dec << std::endl;
      }
      for (auto &needleiter : needle)
      {
        std::cout << "needle: 0x" << std::hex << needleiter << std::dec << std::endl;
      }
    }
  match_again:
    // If found, std::search returns an iterator to the first element of the subsequence
    auto it = std::search(m_HayStack.begin(), m_HayStack.end(), needle.begin(), needle.end());
    if (it != m_HayStack.end())  // haystack and needle have same size - we have a match
    {
      int position = std::distance(m_HayStack.begin(), it);
      if (position > 0)
      {
        if (Verbosity() > 1)
        {
          std::cout << "need to change evt offset of packet " << iter.first << " by "
                    << position << " counts" << std::endl;
        }
        eventoffset[iter.first] = position;
      }
      else
      {
        if (Verbosity() > 1)
        {
          std::cout << "position: " << position << " All good for packet " << iter.first << " with bcodiff " << std::hex << *needle.begin()
                    << " match with " << *m_HayStack.begin() << std::dec << std::endl;
        }
      }
    }
    else
    {
      std::vector needle_copy = needle;
      needle.pop_back();
      // handle 0x0 in needle for first event (if we have an off by one right out of the gate)
      if (*needle.begin() == 0x0)
      {
        needle.erase(needle.begin());
        auto it2 = std::search(m_HayStack.begin(), m_HayStack.end(), needle.begin(), needle.end());
        if (it2 != m_HayStack.end())
        {
          int position = std::distance(m_HayStack.begin(), it2);
          if (Verbosity() > 1)
          {
            std::cout << "Checking for first event, without first element and popped back last one, position: " << position << std::endl;
            for (auto &hayiter : m_HayStack)
            {
              std::cout << "haystack: 0x" << std::hex << hayiter << std::dec << std::endl;
            }
            for (auto &needleiter : needle)
            {
              std::cout << "needle_copy: 0x" << std::hex << needleiter << std::dec << std::endl;
            }
          }
          int shiftby = position - 1;
          if (shiftby > 0)
          {
            eventoffset[iter.first] = shiftby;
          }
        }
      }
      else
      {
        if (needle.size() >= 1)
        {
          //	  needle.pop_back(); already popped back
          // NOLINTNEXTLINE(hicpp-avoid-goto)
          goto match_again;
        }
        // here we have an event which doesn't match
        // first check if the GL1 dropped an event
        // The m_Gl1DroppedEvent set is filled by the GL1 input manager which compares the packet id and the event number
        if (m_Gl1DroppedEvent.find(m_Gl1PacketMap.begin()->first) != m_Gl1DroppedEvent.end())  // dropped gl1 event
        {
          std::cout << "We have a dropped GL1 event" << std::endl;
          m_Gl1DroppedEvent.erase(m_Gl1PacketMap.begin()->first);
          if (!m_MbdPacketMap.empty())
          {
            DropFirstEvent(m_MbdPacketMap);
          }
          if (!m_LL1PacketMap.empty())
          {
            DropFirstEventLL1(m_LL1PacketMap);
          }
          if (!m_HcalPacketMap.empty())
          {
            DropFirstEvent(m_HcalPacketMap);
          }
          if (!m_CemcPacketMap.empty())
          {
            DropFirstEvent(m_CemcPacketMap);
          }
          if (!m_ZdcPacketMap.empty())
          {
            DropFirstEvent(m_ZdcPacketMap);
          }
          if (!m_SEpdPacketMap.empty())
          {
            DropFirstEvent(m_SEpdPacketMap);
          }
          for (auto inputiter : m_NoGl1InputVector)
          {
            inputiter->AdjustEventOffset(-1);
          }
          return -1;
        }
        // what is left is the event after the skip, let's make a crosscheck if the sum of the clockdiff of the previous 2 events does the trick
        if (Verbosity() > 1)
        {
          std::cout << "hay[0]: 0x" << std::hex << m_HayStack[0] << " hay[1]: 0x" << m_HayStack[1] << " sum: 0x" << (m_HayStack[0] + m_HayStack[1]) << std::dec << std::endl;
          std::cout << "needle bdiff: 0x" << std::hex << *needle.begin() << std::dec << std::endl;
        }
        if (*needle.begin() == ((m_HayStack[0] + m_HayStack[1]) & 0xFFFFFFFFU))
        {
          //	  std::cout << "Skipped event" << std::endl;
          eventoffset[iter.first] = 1;
          if (!m_MbdPacketMap.empty())
          {
            AdjustBcoDiff(m_MbdPacketMap, iter.first, m_HayStack[1]);
          }
          if (!m_CemcPacketMap.empty())
          {
            AdjustBcoDiff(m_CemcPacketMap, iter.first, m_HayStack[1]);
          }
          if (!m_HcalPacketMap.empty())
          {
            AdjustBcoDiff(m_HcalPacketMap, iter.first, m_HayStack[1]);
          }

          if (!m_ZdcPacketMap.empty())
          {
            AdjustBcoDiff(m_ZdcPacketMap, iter.first, m_HayStack[1]);
          }
          if (!m_SEpdPacketMap.empty())
          {
            AdjustBcoDiff(m_SEpdPacketMap, iter.first, m_HayStack[1]);
          }
          if (!m_LL1PacketMap.empty())
          {
            AdjustBcoDiffLL1(m_LL1PacketMap, iter.first, m_HayStack[1]);
          }
        }
      }
    }
    //      std::cout << "Sequence found at position: " << position << std::endl;
    if (Verbosity() > 5)
    {
      std::cout << PHWHERE << "haystack size after position check " << m_HayStack.size() << std::endl;
      for (auto iter_1 : m_HayStack)
      {
        std::cout << "haystack: 0x" << std::hex << iter_1 << std::dec << std::endl;
      }
      std::cout << PHWHERE << "needle size of packet " << iter.first << " after position check " << needle.size() << std::endl;
      for (auto &needleiter : needle)
      {
        std::cout << "needle: 0x" << std::hex << needleiter << std::dec << std::endl;
      }
    }
  }

  if (!eventoffset.empty())
  {
    // just loop over all input managers, if it has this packet it will adjust the event number offset for it
    for (auto iter : m_TriggerInputVector)
    {
      for (auto &offiter : eventoffset)
      {
        iter->AdjustEventNumberOffset(offiter.first, offiter.second);
        iter->AdjustPacketMap(offiter.first, offiter.second);
      }
    }
    if (!m_MbdPacketMap.empty())
    {
      ShiftEvents(m_MbdPacketMap, eventoffset, "mbd");
    }
    if (!m_CemcPacketMap.empty())
    {
      ShiftEvents(m_CemcPacketMap, eventoffset, "cemc");
    }
    if (!m_HcalPacketMap.empty())
    {
      ShiftEvents(m_HcalPacketMap, eventoffset, "hcal");
    }

    if (!m_ZdcPacketMap.empty())
    {
      ShiftEvents(m_ZdcPacketMap, eventoffset, "zdc");
    }
    if (!m_SEpdPacketMap.empty())
    {
      ShiftEvents(m_SEpdPacketMap, eventoffset, "zdc");
    }
    if (!m_LL1PacketMap.empty())
    {
      ShiftEventsLL1(m_LL1PacketMap, eventoffset, "ll1");
    }
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::FillNeedle(std::map<int, CaloPacketInfo>::iterator begin, std::map<int, CaloPacketInfo>::iterator end, const std::string &name)
{
  // here we just reset the needle for each packet
  // in principle we only have to remove discarded events but that might be error prone
  // (e.g. what do you do if you have two subsequent events with the same clkdiff)
  // I leave this for later (if at all)
  auto calomapbegin = begin;
  for (auto &pktiter : calomapbegin->second.CaloSinglePacketMap)
  {
    //    std::cout << "Clearing Needle for packet " << pktiter.first << std::endl;
    m_NeedleMap[pktiter.first].clear();
  }
  // this handles the first event where we do not have the bco diff to the previous event
  // for subsequent calls we take the bco diff for the first event in the needle from
  // the cached bco diffs in the bco diffmap
  if (calomapbegin->second.BcoDiffMap.empty())  // This is for the first event, init bco diff to 0x0
  {
    for (auto &pktiter : calomapbegin->second.CaloSinglePacketMap)
    {
      calomapbegin->second.BcoDiffMap[pktiter.first] = 0x0;
      m_NeedleMap[pktiter.first].push_back(0x0);
      // std::cout << "Startup: Pushing 0x0 into packet " << pktiter.first
      // 		<< " for event " << sepdhititer->first << std::endl;
    }
  }
  else
  {
    for (auto &pktiter : calomapbegin->second.CaloSinglePacketMap)
    {
      m_NeedleMap[pktiter.first].push_back(calomapbegin->second.BcoDiffMap[pktiter.first]);
    }
  }
  //      std::cout << PHWHERE << name <<" event: " <<  sepdhititer->first << std::endl;
  for (auto sepdhititer = begin; sepdhititer != end; ++sepdhititer)
  {
    auto nextIt = std::next(sepdhititer);
    if (nextIt != end)
    {
      // This was supposed to save cpu cycles by skipping bcos which have already been set
      // but this needs more thought - if we have event mixing and a resync the BcoDiffMap is partly filled
      // and a simple check if it exists (for any packet) results in further event mixing
      //       if (!nextIt->second.BcoDiffMap.empty())  // this event was already handled, BcoDiffMap is already filled
      //       {
      //         auto pktiter = sepdhititer->second.CaloSinglePacketMap.begin();
      // 	auto checkiter = sepdhititer->second.BcoDiffMap.find(pktiter->first);
      // 	if (checkiter != sepdhititer->second.BcoDiffMap.end())
      // 	{
      // 	  continue;
      // 	}
      // 	std::cout << "could not find bco diff for packet " << pktiter->first << " in event "
      // 		  <<  calomapbegin->first << std::endl;
      // 	// for (auto bcoiter : sepdhititer->second.BcoDiffMap)
      // 	// {
      // 	//   std::cout << "Skipping 0x" << std::hex << bcoiter.second << " into packet " << std::dec << bcoiter.first
      // 	// 	    << " for event " << calomapbegin->first << std::endl;
      // 	// }

      //       }
      std::set<uint64_t> bcodiffs;
      for (auto &pktiter : sepdhititer->second.CaloSinglePacketMap)
      {
        uint64_t prev_bco = pktiter.second->getBCO();
        int prev_packetid = pktiter.first;
        // std::cout << "event " << sepdhititer->first << " packet id: " << prev_packetid
        //  	  << " prev_bco: 0x" << std::hex << prev_bco << std::dec << std::endl;
        auto currpkt = nextIt->second.CaloSinglePacketMap.find(prev_packetid);  //->find(prev_packetid);
        if (currpkt != nextIt->second.CaloSinglePacketMap.end())
        {
          uint64_t curr_bco = currpkt->second->getBCO();
          uint64_t diffbco = curr_bco - prev_bco;
          nextIt->second.BcoDiffMap[prev_packetid] = diffbco;
          if (Verbosity() > 1)
          {
            std::cout << PHWHERE << name << " packet " << prev_packetid << ", prev_bco 0x: " << std::hex
                      << prev_bco << ", curr_bco: 0x" << curr_bco << ", diff: 0x"
                      << diffbco << std::dec << std::endl;
            std::cout << "Pushing 0x" << std::hex << diffbco << " into needle for packet " << std::dec << prev_packetid << std::endl;
          }
          m_NeedleMap[prev_packetid].push_back(diffbco);
          if (Verbosity() > 1)
          {
            std::cout << "Eventnumber " << nextIt->first << ", bco: 0x" << std::hex << pktiter.second->getBCO()
                      << ", diff: 0x" << diffbco << std::dec << std::endl;
          }
          bcodiffs.insert(diffbco);
        }
        else
        {
          if (Verbosity() > 1)
          {
            std::cout << PHWHERE << "Could not find packet " << prev_packetid << " in event " << sepdhititer->first
                      << std::endl;
          }
        }
      }
      if (bcodiffs.size() > 1)
      {
        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << " different bco diffs for " << name << " packets for event " << nextIt->first << std::endl;
        }
      }
    }
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::FillNeedleLL1(std::map<int, LL1PacketInfo>::iterator begin, std::map<int, LL1PacketInfo>::iterator end, const std::string &name)
{
  // here we just reset the needle for each packet
  // in principle we only have to remove discarded events but that might be error prone
  // (e.g. what do you do if you have two subsequent events with the same clkdiff)
  // I leave this for later (if at all)
  auto ll1packetbegin = begin;
  for (auto &pktiter : ll1packetbegin->second.LL1SinglePacketMap)
  {
    m_NeedleMap[pktiter.first].clear();
  }
  // here we refill the needle for every packet with the cached BCO differences
  for (auto sepdhititer = begin; sepdhititer != end; ++sepdhititer)
  {
    for (auto bcoiter : sepdhititer->second.BcoDiffMap)
    {
      m_NeedleMap[bcoiter.first].push_back(bcoiter.second);
    }
  }
  // here we calculate the bco diff to the previous event and update the cached bco difference
  // only for events where we haven't done this yet (check of the bco diff map is empty)
  for (auto sepdhititer = begin; sepdhititer != end; ++sepdhititer)
  {
    if (sepdhititer->second.BcoDiffMap.empty())  // This is for the first event, init bco diff to 0x0
    {
      for (auto &pktiter : sepdhititer->second.LL1SinglePacketMap)
      {
        sepdhititer->second.BcoDiffMap[pktiter.first] = 0x0;
        m_NeedleMap[pktiter.first].push_back(0x0);
        // std::cout << "Startup: Pushing 0x0 into packet " << pktiter.first
        // 		<< " for event " << sepdhititer->first << std::endl;
      }
    }
    auto nextIt = std::next(sepdhititer);
    if (nextIt != end)
    {
      // This was supposed to save cpu cycles by skipping bcos which have already been set
      // but this needs more thought - if we have event mixing and a resync the BcoDiffMap is partly filled
      // and a simple check if it exists (for any packet) results in further event mixing
      // if (!nextIt->second.BcoDiffMap.empty())  // this event was already handled, BcoDiffMap is already filled
      // {
      //   continue;
      // }
      std::set<uint64_t> bcodiffs;
      for (auto &pktiter : sepdhititer->second.LL1SinglePacketMap)
      {
        uint64_t prev_bco = pktiter.second->getBCO();
        int prev_packetid = pktiter.first;
        auto currpkt = nextIt->second.LL1SinglePacketMap.find(prev_packetid);  //->find(prev_packetid);
        if (currpkt != nextIt->second.LL1SinglePacketMap.end())
        {
          uint64_t curr_bco = currpkt->second->getBCO();
          uint64_t diffbco = curr_bco - prev_bco;
          nextIt->second.BcoDiffMap[prev_packetid] = diffbco;
          if (Verbosity() > 11)
          {
            std::cout << PHWHERE << name << " packet " << prev_packetid << ", prev_bco 0x: " << std::hex
                      << prev_bco << ", curr_bco: 0x" << curr_bco << ", diff: 0x"
                      << diffbco << std::dec << std::endl;
            std::cout << "Pushing 0x" << std::hex << diffbco << " into needle for packet " << std::dec << prev_packetid << std::endl;
          }
          m_NeedleMap[prev_packetid].push_back(diffbco);
          bcodiffs.insert(diffbco);
        }
      }
      if (bcodiffs.size() > 1)
      {
        std::cout << PHWHERE << " different bco diffs for " << name << " packets for event " << nextIt->first << std::endl;
      }
    }
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::ShiftEvents(std::map<int, CaloPacketInfo> &PacketInfoMap, std::map<int, int> &eventoffset, const std::string &name)
{
  std::vector<int> eventnumbers;
  std::set<int> packet_ids;
  for (auto pktid : eventoffset)
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "inserting packet " << pktid.first << " to bad boys" << std::endl;
    }
    packet_ids.insert(pktid.first);
  }
  for (auto sepdhititer = PacketInfoMap.rbegin(); sepdhititer != PacketInfoMap.rend(); ++sepdhititer)
  {
    eventnumbers.push_back(sepdhititer->first);
  }
  // we loop over the event numbers instead of the map, since inserting/extracting entries updates the iterators
  // which just breaks the general idea of looping over them once
  for (auto evtnumiter : eventnumbers)
  {
    auto &sepdhititer = PacketInfoMap[evtnumiter];
    for (auto pktiditer : packet_ids)
    {
      //      std::cout << PHWHERE <<  "handling pkt no: " << pktiditer << std::endl;
      auto offsetiter = eventoffset.find(pktiditer);
      if (offsetiter != eventoffset.end())
      {
        int newevent = evtnumiter + offsetiter->second;
        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << name << " moving packet " << pktiditer << " from event " << evtnumiter << " to " << newevent << std::endl;
        }
        auto nh = sepdhititer.CaloSinglePacketMap.extract(pktiditer);
        auto nhbco = sepdhititer.BcoDiffMap.extract(pktiditer);
        //	  std::cout <<  PHWHERE << "size of CaloSinglePacketMap: " <<  PacketInfoMap.size() << std::endl;
        PacketInfoMap[newevent].CaloSinglePacketMap.insert(std::move(nh));
        PacketInfoMap[newevent].BcoDiffMap.insert(std::move(nhbco));
      }
    }
  }
  if (Verbosity() > 1)
  {
    for (auto &sepdeventiter : PacketInfoMap)
    {
      {
        std::cout << "size of map for event " << sepdeventiter.first << " is " << sepdeventiter.second.CaloSinglePacketMap.size() << std::endl;
      }
    }
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::ShiftEventsLL1(std::map<int, LL1PacketInfo> &PacketInfoMap, std::map<int, int> &eventoffset, const std::string &name)
{
  std::vector<int> eventnumbers;
  std::set<int> packet_ids;
  for (auto pktid : eventoffset)
  {
    packet_ids.insert(pktid.first);
  }
  for (auto sepdhititer = PacketInfoMap.rbegin(); sepdhititer != PacketInfoMap.rend(); ++sepdhititer)
  {
    eventnumbers.push_back(sepdhititer->first);
  }
  // we loop over the event numbers instead of the map, since inserting/extracting entries updates the iterators
  // which just breaks the general idea of looping over them once
  for (auto evtnumiter : eventnumbers)
  {
    auto &sepdhititer = PacketInfoMap[evtnumiter];
    //      std::cout << PHWHERE << " Handling event no: " << evtnumiter << " for " << name << std::endl;
    for (auto pktiditer : packet_ids)
    {
      //	std::cout << PHWHERE <<  "handling pkt no: " << pktiditer << std::endl;
      auto offsetiter = eventoffset.find(pktiditer);
      if (offsetiter != eventoffset.end())
      {
        int newevent = evtnumiter + offsetiter->second;
        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << name << " moving packet " << pktiditer << " from event " << evtnumiter << " to " << newevent << std::endl;
        }
        auto nh = sepdhititer.LL1SinglePacketMap.extract(pktiditer);
        auto nhbco = sepdhititer.BcoDiffMap.extract(pktiditer);
        //	  std::cout <<  PHWHERE << "size of LL1SinglePacketMap: " <<  PacketInfoMap.size() << std::endl;
        PacketInfoMap[newevent].LL1SinglePacketMap.insert(std::move(nh));
        PacketInfoMap[newevent].BcoDiffMap.insert(std::move(nhbco));
      }
    }
  }
  if (Verbosity() > 1)
  {
    for (auto &sepdeventiter : PacketInfoMap)
    {
      {
        std::cout << "size of map for event " << sepdeventiter.first << " is " << sepdeventiter.second.LL1SinglePacketMap.size() << std::endl;
      }
    }
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::DropFirstEvent(std::map<int, CaloPacketInfo> &PacketInfoMap)
{
  //  Print("CEMCMAP");
  PacketInfoMap.erase(PacketInfoMap.begin());
  // std::cout << "deleted first event, cemcmap now: " << std::endl;
  // Print("CEMCMAP");

  std::vector<int> events;
  events.reserve(PacketInfoMap.size());
  for (const auto &packetloop : PacketInfoMap)
  {
    events.push_back(packetloop.first);
  }
  for (auto evtiter : events)
  {
    //    std::cout << "moving event " << evtiter << " to " << (evtiter - 1) << std::endl;
    auto nh = PacketInfoMap[evtiter];
    //    PacketInfoMap.insert(std::make_pair(evtiter-1,std::move(nh)));
    PacketInfoMap[evtiter - 1] = std::move(nh);
    PacketInfoMap.erase(evtiter);
    //    Print("CEMCMAP");
  }
  // std::cout << "before killing last entry" << std::endl;
  // Print("CEMCMAP");
  // std::cout << "That is it" << std::endl;
  // Print("CEMCMAP");
  return 0;
}

int Fun4AllPrdfInputTriggerManager::DropFirstEventLL1(std::map<int, LL1PacketInfo> &PacketInfoMap)
{
  //  Print("CEMCMAP");
  PacketInfoMap.erase(PacketInfoMap.begin());
  // std::cout << "deleted first event, cemcmap now: " << std::endl;
  // Print("CEMCMAP");

  std::vector<int> events;
  events.reserve(PacketInfoMap.size());
  for (const auto &packetloop : PacketInfoMap)
  {
    events.push_back(packetloop.first);
  }
  for (auto evtiter : events)
  {
    //    std::cout << "moving event " << evtiter << " to " << (evtiter - 1) << std::endl;
    auto nh = PacketInfoMap[evtiter];
    //    PacketInfoMap.insert(std::make_pair(evtiter-1,std::move(nh)));
    PacketInfoMap[evtiter - 1] = std::move(nh);
    PacketInfoMap.erase(evtiter);
    //    Print("CEMCMAP");
  }
  // std::cout << "before killing last entry" << std::endl;
  // Print("CEMCMAP");
  // std::cout << "That is it" << std::endl;
  // Print("CEMCMAP");
  return 0;
}

int Fun4AllPrdfInputTriggerManager::AdjustBcoDiff(std::map<int, CaloPacketInfo> &PacketInfoMap, int packetid, uint64_t bcodiff)
{
  auto calomapiter = PacketInfoMap.begin();
  auto pkt = calomapiter->second.BcoDiffMap.find(packetid);
  if (pkt != calomapiter->second.BcoDiffMap.end())
  {
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " changing BCO for packet " << packetid << " for event " << calomapiter->first
                << " to 0x" << std::hex << bcodiff << std::dec << std::endl;
    }
    pkt->second = bcodiff;
  }
  return 0;
}

int Fun4AllPrdfInputTriggerManager::AdjustBcoDiffLL1(std::map<int, LL1PacketInfo> &PacketInfoMap, int packetid, uint64_t bcodiff)
{
  auto calomapiter = PacketInfoMap.begin();
  auto pkt = calomapiter->second.BcoDiffMap.find(packetid);
  if (pkt != calomapiter->second.BcoDiffMap.end())
  {
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " changing BCO for packet " << packetid << " for event " << calomapiter->first
                << " to 0x" << std::hex << bcodiff << std::dec << std::endl;
    }
    pkt->second = bcodiff;
  }
  return 0;
}

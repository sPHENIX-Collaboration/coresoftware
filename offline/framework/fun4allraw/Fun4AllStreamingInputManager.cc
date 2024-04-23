#include "Fun4AllStreamingInputManager.h"

#include "InputManagerType.h"
#include "SingleStreamingInput.h"

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>
#include <ffarawobjects/MvtxRawEvtHeader.h>
#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <ffaobjects/SyncObject.h>  // for SyncObject
#include <ffaobjects/SyncObjectv1.h>

#include <frog/FROG.h>

#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <algorithm>  // for max
#include <cassert>
#include <cstdint>  // for uint64_t, uint16_t
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllStreamingInputManager::Fun4AllStreamingInputManager(const std::string &name, const std::string &dstnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, dstnodename, topnodename)
  , m_SyncObject(new SyncObjectv1())
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());
  return;
}

Fun4AllStreamingInputManager::~Fun4AllStreamingInputManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete m_SyncObject;
  // clear leftover raw event maps and vectors with poolreaders
  // GL1
  for (auto iter : m_Gl1InputVector)
  {
    delete iter;
  }

  m_Gl1InputVector.clear();

  // MVTX
  for (auto const &mapiter : m_MvtxRawHitMap)
  {
    for (auto mvtxhititer : mapiter.second.MvtxRawHitVector)
    {
      delete mvtxhititer;
    }
  }
  m_MvtxRawHitMap.clear();

  for (auto iter : m_MvtxInputVector)
  {
    delete iter;
  }
  m_MvtxInputVector.clear();

  // INTT
  for (auto const &mapiter : m_InttRawHitMap)
  {
    for (auto intthititer : mapiter.second.InttRawHitVector)
    {
      delete intthititer;
    }
  }
  m_InttRawHitMap.clear();

  for (auto iter : m_InttInputVector)
  {
    delete iter;
  }
  m_InttInputVector.clear();

  // TPC
  for (auto const &mapiter : m_TpcRawHitMap)
  {
    for (auto tpchititer : mapiter.second.TpcRawHitVector)
    {
      delete tpchititer;
    }
  }

  m_TpcRawHitMap.clear();
  for (auto iter : m_TpcInputVector)
  {
    delete iter;
  }
  m_TpcInputVector.clear();

  // Micromegas

  for (auto const &mapiter : m_MicromegasRawHitMap)
  {
    for (auto micromegashititer : mapiter.second.MicromegasRawHitVector)
    {
      delete micromegashititer;
    }
  }
  for (auto iter : m_MicromegasInputVector)
  {
    delete iter;
  }

  m_MicromegasInputVector.clear();
}

int Fun4AllStreamingInputManager::run(const int /*nevents*/)
{
  int iret = 0;
  if (m_gl1_registered_flag)  // Gl1 first to get the reference
  {
    iret += FillGl1();
  }
  if (m_intt_registered_flag)
  {
    iret += FillIntt();
  }
  if (m_mvtx_registered_flag)
  {
    iret += FillMvtx();
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

int Fun4AllStreamingInputManager::fileclose()
{
  return 0;
}

void Fun4AllStreamingInputManager::Print(const std::string &what) const
{
  if (what == "TPC")
  {
    for (auto &iter : m_TpcRawHitMap)
    {
      std::cout << "bco: " << std::hex << iter.first << std::dec << std::endl;
      for (auto &itervec : iter.second.TpcRawHitVector)
      {
        std::cout << "hit: " << std::hex << itervec << std::dec << std::endl;
        itervec->identify();
      }
    }
  }
  if (what == "ALL" || what == "INPUTFILES")
  {
    std::cout << "-----------------------------" << std::endl;
    for (const auto &iter : m_Gl1InputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
    for (const auto &iter : m_MvtxInputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
    for (const auto &iter : m_InttInputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
    for (const auto &iter : m_TpcInputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
    for (const auto &iter : m_MicromegasInputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
  }
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllStreamingInputManager::ResetEvent()
{
  // for (auto iter : m_EvtInputVector)
  // {
  //   iter->CleanupUsedPackets(m_CurrentBeamClock);
  // }
  //  m_SyncObject->Reset();
  return 0;
}

int Fun4AllStreamingInputManager::PushBackEvents(const int /*i*/)
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
  //        << " Fun4AllStreamingInputManager cannot push back " << i << " events into file"
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

int Fun4AllStreamingInputManager::GetSyncObject(SyncObject **mastersync)
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

int Fun4AllStreamingInputManager::SyncIt(const SyncObject *mastersync)
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

std::string Fun4AllStreamingInputManager::GetString(const std::string &what) const
{
  std::cout << PHWHERE << " called with " << what << " , returning empty string" << std::endl;
  return "";
}

void Fun4AllStreamingInputManager::registerStreamingInput(SingleStreamingInput *evtin, InputManagerType::enu_subsystem system)
{
  evtin->StreamingInputManager(this);
  evtin->CreateDSTNode(m_topNode);
  evtin->ConfigureStreamingInputManager();
  switch (system)
  {
  case InputManagerType::MVTX:
    m_mvtx_registered_flag = true;
    m_MvtxInputVector.push_back(evtin);
    break;
  case InputManagerType::INTT:
    m_intt_registered_flag = true;
    m_InttInputVector.push_back(evtin);
    break;
  case InputManagerType::TPC:
    m_tpc_registered_flag = true;
    m_TpcInputVector.push_back(evtin);
    break;
  case InputManagerType::MICROMEGAS:
    m_micromegas_registered_flag = true;
    m_MicromegasInputVector.push_back(evtin);
    break;
  case InputManagerType::GL1:
    m_gl1_registered_flag = true;
    m_Gl1InputVector.push_back(evtin);
    break;
  default:
    std::cout << "invalid subsystem flag " << system << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  if (Verbosity() > 3)
  {
    std::cout << "registering " << evtin->Name()
              << " number of registered inputs: "
              << m_Gl1InputVector.size() + m_InttInputVector.size() + m_MicromegasInputVector.size() + m_MvtxInputVector.size() + m_TpcInputVector.size()
              << std::endl;
  }
}

void Fun4AllStreamingInputManager::AddGl1RawHit(uint64_t bclk, Gl1RawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding gl1 hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_Gl1RawHitMap[bclk].Gl1RawHitVector.push_back(hit);
}

void Fun4AllStreamingInputManager::AddMvtxRawHit(uint64_t bclk, MvtxRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding mvtx hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_MvtxRawHitMap[bclk].MvtxRawHitVector.push_back(hit);
}

void Fun4AllStreamingInputManager::AddMvtxFeeId(uint64_t bclk, uint16_t feeid)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding mvtx feeid to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_MvtxRawHitMap[bclk].MvtxFeeIds.insert(feeid);
}

void Fun4AllStreamingInputManager::AddMvtxL1TrgBco(uint64_t bclk, uint64_t lv1Bco)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding mvtx L1Trg to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_MvtxRawHitMap[bclk].MvtxL1TrgBco.insert(lv1Bco);
}

void Fun4AllStreamingInputManager::AddInttRawHit(uint64_t bclk, InttRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding intt hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_InttRawHitMap[bclk].InttRawHitVector.push_back(hit);
  m_InttPacketFeeBcoMap[hit->get_packetid()][hit->get_fee()] = bclk;
}

void Fun4AllStreamingInputManager::AddMicromegasRawHit(uint64_t bclk, MicromegasRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding micromegas hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_MicromegasRawHitMap[bclk].MicromegasRawHitVector.push_back(hit);
}

void Fun4AllStreamingInputManager::AddTpcRawHit(uint64_t bclk, TpcRawHit *hit)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding tpc hit to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  m_TpcRawHitMap[bclk].TpcRawHitVector.push_back(hit);
}

int Fun4AllStreamingInputManager::FillGl1()
{
  // unsigned int alldone = 0;
  for (auto iter : m_Gl1InputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillGl1 - fill pool for " << iter->Name() << std::endl;
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
  if (m_Gl1RawHitMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  //    std::cout << "stashed gl1 BCOs: " << m_Gl1RawHitMap.size() << std::endl;
  Gl1RawHit *gl1rawhit = findNode::getClass<Gl1RawHit>(m_topNode, "GL1RAWHIT");
  //  std::cout << "before filling m_Gl1RawHitMap size: " <<  m_Gl1RawHitMap.size() << std::endl;
  for (auto gl1hititer : m_Gl1RawHitMap.begin()->second.Gl1RawHitVector)
  {
    if (Verbosity() > 1)
    {
      gl1hititer->identify();
    }
    gl1rawhit->CopyFrom(gl1hititer);
    MySyncManager()->CurrentEvent(gl1rawhit->getEvtSequence());
    m_RefBCO = gl1hititer->get_bco();
    m_RefBCO = m_RefBCO & 0xFFFFFFFFFFU;  // 40 bits (need to handle rollovers)
  }
  for (auto iter : m_Gl1InputVector)
  {
    iter->CleanupUsedPackets(m_Gl1RawHitMap.begin()->first);
  }
  m_Gl1RawHitMap.begin()->second.Gl1RawHitVector.clear();
  m_Gl1RawHitMap.erase(m_Gl1RawHitMap.begin());
  // std::cout << "size  m_Gl1RawHitMap: " <<  m_Gl1RawHitMap.size()
  // 	    << std::endl;
  return 0;
}

int Fun4AllStreamingInputManager::FillIntt()
{
  int iret = FillInttPool();
  if (iret)
  {
    return iret;
  }
  // unsigned int alldone = 0;
  //     std::cout << "stashed intt BCOs: " << m_InttRawHitMap.size() << std::endl;
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(m_topNode, "INTTRAWHIT");
  //  std::cout << "before filling m_InttRawHitMap size: " <<  m_InttRawHitMap.size() << std::endl;
  // !m_InttRawHitMap.empty() is implicitely handled and the check is expensive
  // FillInttPool() contains this check already and will return non zero
  // so here m_InttRawHitMap will always contain entries
  uint64_t select_crossings = m_intt_bco_range;
  if (m_RefBCO > 0)
  {
    select_crossings += m_RefBCO;
  }
  else
  {
    select_crossings += m_InttRawHitMap.begin()->first;
  }
  if (Verbosity() > 2)
  {
    std::cout << "select INTT crossings"
              << " from 0x" << std::hex << m_RefBCO - m_intt_negative_bco
              << " to 0x" << select_crossings - m_intt_negative_bco
              << " for ref BCO " << m_RefBCO
              << std::dec << std::endl;
  }
  while (m_InttRawHitMap.begin()->first < m_RefBCO - m_intt_negative_bco)
  {
    std::cout << "Intt BCO: 0x" << std::hex << m_InttRawHitMap.begin()->first
              << " corrected for negative offset: 0x" << m_InttRawHitMap.begin()->first + m_intt_negative_bco
              << " smaller than GL1 BCO: 0x" << m_RefBCO
              << " corrected for range: 0x" << select_crossings
              << ", ditching this bco" << std::dec << std::endl;
    for (auto iter : m_InttInputVector)
    {
      iter->CleanupUsedPackets(m_InttRawHitMap.begin()->first);
    }
    m_InttRawHitMap.begin()->second.InttRawHitVector.clear();
    m_InttRawHitMap.erase(m_InttRawHitMap.begin());
    iret = FillInttPool();
    if (iret)
    {
      return iret;
    }
  }
  while (m_InttRawHitMap.begin()->first <= select_crossings - m_intt_negative_bco)
  {
    for (auto intthititer : m_InttRawHitMap.begin()->second.InttRawHitVector)
    {
      if (Verbosity() > 1)
      {
        intthititer->identify();
      }
      inttcont->AddHit(intthititer);
    }
    for (auto iter : m_InttInputVector)
    {
      iter->CleanupUsedPackets(m_InttRawHitMap.begin()->first);
    }
    m_InttRawHitMap.begin()->second.InttRawHitVector.clear();
    m_InttRawHitMap.erase(m_InttRawHitMap.begin());
    if (m_InttRawHitMap.empty())
    {
      break;
    }
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillMvtx()
{
  int iret = FillMvtxPool();
  if (iret)
  {
    return iret;
  }
  MvtxRawEvtHeader *mvtxEvtHeader = findNode::getClass<MvtxRawEvtHeader>(m_topNode, "MVTXRAWEVTHEADER");
  if (!mvtxEvtHeader)
  {
    std::cout << "ERROR: MVTXRAWEVTHEADER node not found, exit. " << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  MvtxRawHitContainer *mvtxcont = findNode::getClass<MvtxRawHitContainer>(m_topNode, "MVTXRAWHIT");
  if (!mvtxcont)
  {
    std::cout << "ERROR: MVTXRAWHIT node not found, exit. " << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  // std::cout << "before filling m_MvtxRawHitMap size: " <<  m_MvtxRawHitMap.size() << std::endl;
  uint64_t select_crossings = m_mvtx_bco_range;
  if (m_RefBCO > 0)
  {
    select_crossings += m_RefBCO;
  }
  else
  {
    select_crossings += m_MvtxRawHitMap.begin()->first;
  }

  uint64_t ref_bco_minus_range = m_RefBCO < m_mvtx_bco_range ? 0 : m_RefBCO - m_mvtx_bco_range;
  if (Verbosity() > 2)
  {
    std::cout << "select MVTX crossings"
              << " from 0x" << std::hex << ref_bco_minus_range
              << " to 0x" << select_crossings - m_mvtx_bco_range
              << std::dec << std::endl;
  }
  // m_MvtxRawHitMap.empty() does not need to be checked here, FillMvtxPool returns non zero
  // if this map is empty which is handled above
  //All three values used in the while loop evaluation are unsigned ints. If m_RefBCO is < m_mvtx_bco_range then we will overflow and delete all hits
  while (m_MvtxRawHitMap.begin()->first < ref_bco_minus_range)
  {
    if (Verbosity() > 2)
    {
      std::cout << "ditching mvtx bco 0x" << std::hex << m_MvtxRawHitMap.begin()->first << ", ref: 0x" << m_RefBCO << std::dec << std::endl;
    }
    for (auto iter : m_MvtxInputVector)
    {
      iter->CleanupUsedPackets(m_MvtxRawHitMap.begin()->first);
    }
    m_MvtxRawHitMap.begin()->second.MvtxRawHitVector.clear();
    m_MvtxRawHitMap.erase(m_MvtxRawHitMap.begin());
    iret = FillMvtxPool();
    if (iret)
    {
      return iret;
    }
  }
  // again m_MvtxRawHitMap.empty() is handled by return of FillMvtxPool()
  if (Verbosity() > 2)
  {
    std::cout << "after ditching, mvtx bco: 0x" << std::hex << m_MvtxRawHitMap.begin()->first << ", ref: 0x" << m_RefBCO
              << std::dec << std::endl;
  }
  auto mvtxRawHitInfoIt = m_MvtxRawHitMap.begin();

  mvtxEvtHeader->AddFeeId(mvtxRawHitInfoIt->second.MvtxFeeIds);
  mvtxEvtHeader->AddL1Trg(mvtxRawHitInfoIt->second.MvtxL1TrgBco);
  while (m_MvtxRawHitMap.begin()->first <= select_crossings - m_mvtx_bco_range)
  {
    if (Verbosity() > 2)
    {
      std::cout << "Adding 0x" << std::hex << m_MvtxRawHitMap.begin()->first
                << " ref: 0x" << select_crossings << std::dec << std::endl;
    }
    for (auto mvtxhititer : m_MvtxRawHitMap.begin()->second.MvtxRawHitVector)
    {
      if (Verbosity() > 1)
      {
        mvtxhititer->identify();
      }
      mvtxcont->AddHit(mvtxhititer);
    }
    for (auto iter : m_MvtxInputVector)
    {
      iter->CleanupUsedPackets(m_MvtxRawHitMap.begin()->first);
    }
    m_MvtxRawHitMap.begin()->second.MvtxRawHitVector.clear();
    m_MvtxRawHitMap.erase(m_MvtxRawHitMap.begin());
    // m_MvtxRawHitMap.empty() need to be checked here since we do not call FillPoolMvtx()
    if (m_MvtxRawHitMap.empty())
    {
      break;
    }
  }

  return 0;
}

//_______________________________________________________
int Fun4AllStreamingInputManager::FillMicromegas()
{
  int iret = FillMicromegasPool();
  if (iret)
  {
    return iret;
  }

  auto container = findNode::getClass<MicromegasRawHitContainer>(m_topNode, "MICROMEGASRAWHIT");
  uint64_t select_crossings = m_micromegas_bco_range;
  if (m_RefBCO > 0)
  {
    select_crossings += m_RefBCO;
  }
  else
  {
    select_crossings += m_MicromegasRawHitMap.begin()->first;
  }
  if (Verbosity() > 2)
  {
    std::cout << "select MicroMegas crossings"
              << " from 0x" << std::hex << m_RefBCO - m_micromegas_bco_range
              << " to 0x" << select_crossings - m_micromegas_bco_range
              << " for ref BCO " << m_RefBCO
              << std::dec << std::endl;
  }

  // m_MicromegasRawHitMap.empty() does not need to be checked here, FillMicromegasPool returns non zero
  // if this map is empty which is handled above
  while (m_MicromegasRawHitMap.begin()->first < m_RefBCO - m_micromegas_negative_bco)
  {
    std::cout << "Micromegas BCO: 0x" << std::hex << m_MicromegasRawHitMap.begin()->first
              << " smaller than GL1 BCO: 0x" << m_RefBCO
              << ", ditching this bco" << std::dec << std::endl;
    for (auto iter : m_MicromegasInputVector)
    {
      iter->CleanupUsedPackets(m_MicromegasRawHitMap.begin()->first);
    }
    m_MicromegasRawHitMap.begin()->second.MicromegasRawHitVector.clear();
    m_MicromegasRawHitMap.erase(m_MicromegasRawHitMap.begin());
    iret = FillMicromegasPool();
    if (iret)
    {
      return iret;
    }
  }

  while ((m_MicromegasRawHitMap.begin()->first) <= select_crossings - m_micromegas_negative_bco)
  {
    for (const auto &hititer : m_MicromegasRawHitMap.begin()->second.MicromegasRawHitVector)
    {
      container->AddHit(hititer);
    }

    for (const auto &iter : m_MicromegasInputVector)
    {
      iter->CleanupUsedPackets(m_MicromegasRawHitMap.begin()->first);
    }

    m_MicromegasRawHitMap.begin()->second.MicromegasRawHitVector.clear();
    m_MicromegasRawHitMap.erase(m_MicromegasRawHitMap.begin());
    if (m_MicromegasRawHitMap.empty())
    {
      break;
    }
  }

  return 0;
}

int Fun4AllStreamingInputManager::FillTpc()
{
  int iret = FillTpcPool();
  if (iret)
  {
    return iret;
  }
  TpcRawHitContainer *tpccont = findNode::getClass<TpcRawHitContainer>(m_topNode, "TPCRAWHIT");
  //  std::cout << "before filling m_TpcRawHitMap size: " <<  m_TpcRawHitMap.size() << std::endl;
  uint64_t select_crossings = m_tpc_bco_range;
  if (m_RefBCO > 0)
  {
    select_crossings += m_RefBCO;
  }
  else
  {
    select_crossings += m_TpcRawHitMap.begin()->first;
  }
  if (Verbosity() > 2)
  {
    std::cout << "select TPC crossings"
              << " from 0x" << std::hex << m_RefBCO - m_tpc_negative_bco
              << " to 0x" << select_crossings - m_tpc_negative_bco
              << std::dec << std::endl;
  }
  // m_TpcRawHitMap.empty() does not need to be checked here, FillTpcPool returns non zero
  // if this map is empty which is handled above
  while (m_TpcRawHitMap.begin()->first < m_RefBCO - m_tpc_negative_bco)
  {
    for (auto iter : m_TpcInputVector)
    {
      iter->CleanupUsedPackets(m_TpcRawHitMap.begin()->first);
    }
    m_TpcRawHitMap.begin()->second.TpcRawHitVector.clear();
    m_TpcRawHitMap.erase(m_TpcRawHitMap.begin());
    iret = FillTpcPool();
    if (iret)
    {
      return iret;
    }
  }
  // again m_TpcRawHitMap.empty() is handled by return of FillTpcPool()
  while (m_TpcRawHitMap.begin()->first <= select_crossings - m_tpc_negative_bco)
  {
    for (auto tpchititer : m_TpcRawHitMap.begin()->second.TpcRawHitVector)
    {
      if (Verbosity() > 1)
      {
        tpchititer->identify();
      }
      tpccont->AddHit(tpchititer);
    }
    for (auto iter : m_TpcInputVector)
    {
      iter->CleanupUsedPackets(m_TpcRawHitMap.begin()->first);
    }
    m_TpcRawHitMap.begin()->second.TpcRawHitVector.clear();
    m_TpcRawHitMap.erase(m_TpcRawHitMap.begin());
    if (m_TpcRawHitMap.empty())
    {
      break;
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << "tpc container size: " << tpccont->get_nhits();
    std::cout << ", size  m_TpcRawHitMap: " << m_TpcRawHitMap.size()
              << std::endl;
  }
  if (tpccont->get_nhits() > 500000)
  {
    std::cout << "Resetting TPC Container with number of entries " << tpccont->get_nhits() << std::endl;
    tpccont->Reset();
  }
  return 0;
}

void Fun4AllStreamingInputManager::SetInttBcoRange(const unsigned int i)
{
  m_intt_bco_range = std::max(i, m_intt_bco_range);
}

void Fun4AllStreamingInputManager::SetInttNegativeBco(const unsigned int i)
{
  m_intt_negative_bco = std::max(i, m_intt_negative_bco);
}

void Fun4AllStreamingInputManager::SetMicromegasBcoRange(const unsigned int i)
{
  m_micromegas_bco_range = std::max(i, m_micromegas_bco_range);
}

void Fun4AllStreamingInputManager::SetMicromegasNegativeBco(const unsigned int i)
{
  m_micromegas_negative_bco = std::max(i, m_micromegas_negative_bco);
}

void Fun4AllStreamingInputManager::SetMvtxNegativeBco(const unsigned int i)
{
  m_mvtx_negative_bco = std::max(i, m_mvtx_negative_bco);
}

void Fun4AllStreamingInputManager::SetTpcBcoRange(const unsigned int i)
{
  m_tpc_bco_range = std::max(i, m_tpc_bco_range);
}

void Fun4AllStreamingInputManager::SetTpcNegativeBco(const unsigned int i)
{
  m_tpc_negative_bco = std::max(i, m_tpc_negative_bco);
}

void Fun4AllStreamingInputManager::SetMvtxBcoRange(const unsigned int i)
{
  m_mvtx_bco_range = std::max(i, m_mvtx_bco_range);
}

int Fun4AllStreamingInputManager::FillInttPool()
{
  for (auto iter : m_InttInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillInttPool - fill pool for " << iter->Name() << std::endl;
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
  if (m_InttRawHitMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillTpcPool()
{
  for (auto iter : m_TpcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillTpcPool - fill pool for " << iter->Name() << std::endl;
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
  if (m_TpcRawHitMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillMicromegasPool()
{
  for (auto iter : m_MicromegasInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillMicromegasPool - fill pool for " << iter->Name() << std::endl;
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
  if (m_MicromegasRawHitMap.empty())
  {
    std::cout << "Micromegas are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillMvtxPool()
{
  for (auto iter : m_MvtxInputVector)
  {
    if (Verbosity() > 3)
    {
      std::cout << "Fun4AllStreamingInputManager::FillMvtxPool - fill pool for " << iter->Name() << std::endl;
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
  if (m_MvtxRawHitMap.empty())
  {
    std::cout << "we are done" << std::endl;
    return -1;
  }
  return 0;
}

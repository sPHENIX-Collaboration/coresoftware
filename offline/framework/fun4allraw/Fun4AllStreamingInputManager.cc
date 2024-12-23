#include "Fun4AllStreamingInputManager.h"

#include "InputManagerType.h"
#include "MvtxRawDefs.h"
#include "SingleMicromegasPoolInput.h"
#include "SingleMicromegasPoolInput_v2.h"
#include "SingleMvtxPoolInput.h"
#include "SingleStreamingInput.h"

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>
#include <ffarawobjects/MvtxFeeIdInfov1.h>
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

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <frog/FROG.h>

#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <boost/format.hpp>

#include <TH1.h>
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

  createQAHistos();

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
    for (auto mvtxFeeIdInfo : mapiter.second.MvtxFeeIdInfoVector)
    {
      delete mvtxFeeIdInfo;
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
  m_RefBCO = 0;
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
  // if the streaming flag is set, we only want the first event from the GL1 to
  // get the starting BCO of that run which enables us to dump all the junk which
  // is taken before the run starts in the streaming systems. But we don't want the
  // GL1 in the output, so we do not create its dst node if running in streaming
  if (system == InputManagerType::GL1)
  {
    if (!m_StreamingFlag)
    {
      evtin->CreateDSTNode(m_topNode);
    }
  }
  else
  {
    evtin->CreateDSTNode(m_topNode);
  }
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
    evtin->createQAHistos();
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

void Fun4AllStreamingInputManager::AddGl1RawHit(uint64_t bclk, Gl1Packet *hit)
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

void Fun4AllStreamingInputManager::AddMvtxFeeIdInfo(uint64_t bclk, uint16_t feeid, uint32_t detField)
{
  if (Verbosity() > 1)
  {
    std::cout << "Adding mvtx feeid info to bclk 0x"
              << std::hex << bclk << std::dec << std::endl;
  }
  MvtxFeeIdInfo *feeidInfo = new MvtxFeeIdInfov1();
  feeidInfo->set_bco(bclk);
  feeidInfo->set_feeId(feeid);
  feeidInfo->set_detField(detField);
  m_MvtxRawHitMap[bclk].MvtxFeeIdInfoVector.push_back(feeidInfo);
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
    std::cout << "Gl1RawHitMap is empty - we are done" << std::endl;
    return -1;
  }
  //    std::cout << "stashed gl1 BCOs: " << m_Gl1RawHitMap.size() << std::endl;
  Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(m_topNode, "GL1RAWHIT");
  //  std::cout << "before filling m_Gl1RawHitMap size: " <<  m_Gl1RawHitMap.size() << std::endl;
  for (auto gl1hititer : m_Gl1RawHitMap.begin()->second.Gl1RawHitVector)
  {
    if (Verbosity() > 1)
    {
      gl1hititer->identify();
    }
    if (!m_StreamingFlag)  // if streaming flag is set, the gl1packet is a nullptr
    {
      gl1packet->FillFrom(gl1hititer);
      MySyncManager()->CurrentEvent(gl1packet->getEvtSequence());
    }
    m_RefBCO = gl1hititer->getBCO();
    m_RefBCO = m_RefBCO & 0xFFFFFFFFFFU;  // 40 bits (need to handle rollovers)
                                          //    std::cout << "BCOis " << std::hex << m_RefBCO << std::dec << std::endl;
  }
  // if we run streaming, we only need the first gl1 bco to skip over all the junk
  // which is taken before the daq actually starts. But once we have the first event
  // and set the refBCO to the beginning of the run, we don't want the gl1 anymore
  // so we delete its input manager(s) and unregister it
  // deleting it also deletes all its allocated memory, so we don't have to worry
  // about clearing all gl1 related maps
  if (m_StreamingFlag)
  {
    for (auto iter : m_Gl1InputVector)
    {
      delete iter;
    }
      m_gl1_registered_flag = false;
      m_Gl1InputVector.clear();
  }
  else
  {
    for (auto iter : m_Gl1InputVector)
    {
      iter->CleanupUsedPackets(m_Gl1RawHitMap.begin()->first);
    }
    m_Gl1RawHitMap.begin()->second.Gl1RawHitVector.clear();
    m_Gl1RawHitMap.erase(m_Gl1RawHitMap.begin());
  }
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
  if (!inttcont)
  {
    inttcont = findNode::getClass<InttRawHitContainer>(m_topNode, (*(m_InttInputVector.begin()))->getHitContainerName());
    if (!inttcont)
    {
      std::cout << PHWHERE << "Could not find InttRawHitContainer node in topNode" << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }
  //  std::cout << "before filling m_InttRawHitMap size: " <<  m_InttRawHitMap.size() << std::endl;
  // !m_InttRawHitMap.empty() is implicitely handled and the check is expensive
  // FillInttPool() contains this check already and will return non zero
  // so here m_InttRawHitMap will always contain entries
  uint64_t select_crossings = m_intt_bco_range;
  if (m_RefBCO == 0)
  {
    m_RefBCO = m_InttRawHitMap.begin()->first;
  }
  select_crossings += m_RefBCO;
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
    if (Verbosity() > 2)
    {
      std::cout << "Intt BCO: 0x" << std::hex << m_InttRawHitMap.begin()->first
                << " corrected for negative offset: 0x" << m_InttRawHitMap.begin()->first + m_intt_negative_bco
                << " smaller than GL1 BCO: 0x" << m_RefBCO
                << " corrected for range: 0x" << select_crossings
                << std::dec << " diff: " << (m_RefBCO - m_InttRawHitMap.begin()->first)
                << ", ditching this bco" << std::dec << std::endl;
    }
    for (auto iter : m_InttInputVector)
    {
      iter->CleanupUsedPackets(m_InttRawHitMap.begin()->first);
      iter->clearPacketBClkStackMap(m_InttRawHitMap.begin()->first);
      iter->clearFeeGTML1BCOMap(m_InttRawHitMap.begin()->first);
    }
   
    m_InttRawHitMap.begin()->second.InttRawHitVector.clear();
    m_InttRawHitMap.erase(m_InttRawHitMap.begin());
    
    iret = FillInttPool();
    if (iret)
    {
      return iret;
    }
  }

  unsigned int refbcobitshift = m_RefBCO & 0x3FU;
  h_refbco_intt->Fill(refbcobitshift);
  bool allpackets = true;
  int allpacketsallfees = 0;
  for (auto &p : m_InttInputVector)
  {
    
    // this is on a per packet basis
    auto bcl_stack = p->BclkStackMap();
    auto feebclstack = p->getFeeGTML1BCOMap();
    int packet_id = bcl_stack.begin()->first;
    int histo_to_fill = (packet_id % 10) - 1;
  
    std::set<int> feeidset;
    int fee = 0;
    for (auto &[feeid, gtmbcoset] : feebclstack)
    {
      for (auto &bcl : gtmbcoset)
      {
        auto diff = (m_RefBCO > bcl) ? m_RefBCO - bcl : bcl - m_RefBCO;
        if (diff < 120) { // diff is 1 strobe length of 120 crossings
          h_gl1taggedfee_intt[histo_to_fill][fee]->Fill(refbcobitshift);
          feeidset.insert(feeid);
          // this fee was tagged, go to the next one
          break;
        }
      }
      fee++;
    }
    
    if (feeidset.size() == 14)
    {
      allpacketsallfees++;
      h_taggedAllFees_intt[histo_to_fill]->Fill(refbcobitshift);
    }
    feeidset.clear();
    bool thispacket = false;
        
    for (auto &[packetid, gtmbcoset] : bcl_stack)
    {
      for (auto &gtmbco : gtmbcoset)
      {
        auto diff = (m_RefBCO > gtmbco) ? m_RefBCO - gtmbco : gtmbco - m_RefBCO;
        if (diff < 120)  //diff is 1 strobe length of 120 crossings
        {
          thispacket = true;
          h_gl1tagged_intt[histo_to_fill]->Fill(refbcobitshift);
          break;
        }
      }
    }
    
    if (thispacket == false)
    {
      allpackets = false;
    }
  }
  if (allpackets)
  {
    h_taggedAll_intt->Fill(refbcobitshift);
  }
  if (allpacketsallfees == (int) m_InttInputVector.size())
  {
    h_taggedAllFee_intt->Fill(refbcobitshift);
  }

  //  std::cout << "Checking diff " << (m_InttRawHitMap.begin()->first - (select_crossings - m_intt_negative_bco)) << std::endl;
  while (m_InttRawHitMap.begin()->first <= select_crossings - m_intt_negative_bco)
  {
    for (auto intthititer : m_InttRawHitMap.begin()->second.InttRawHitVector)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Adding intt hit with bco 0x" << std::hex
                  << intthititer->get_bco() << std::dec << std::endl;
        //        intthititer->identify();
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
    mvtxEvtHeader = findNode::getClass<MvtxRawEvtHeader>(m_topNode, (static_cast<SingleMvtxPoolInput *>(*(m_MvtxInputVector.begin())))->getRawEventHeaderName());
    if (!mvtxEvtHeader)
    {
      std::cout << PHWHERE << "ERROR: MVTXRAWEVTHEADER node not found, exit. " << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  MvtxRawHitContainer *mvtxcont = findNode::getClass<MvtxRawHitContainer>(m_topNode, "MVTXRAWHIT");
  if (!mvtxcont)
  {
    mvtxcont = findNode::getClass<MvtxRawHitContainer>(m_topNode, (*(m_MvtxInputVector.begin()))->getHitContainerName());
    if (!mvtxcont)
    {
      std::cout << PHWHERE << "ERROR: MVTXRAWHIT node not found, exit. " << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }
  // std::cout << "before filling m_MvtxRawHitMap size: " <<  m_MvtxRawHitMap.size() << std::endl;
  uint64_t select_crossings = m_mvtx_is_triggered ? 0 : m_mvtx_bco_range;
  if (m_RefBCO == 0)
  {
    m_RefBCO = m_MvtxRawHitMap.begin()->first;
  }
  select_crossings += m_RefBCO;

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
  // All three values used in the while loop evaluation are unsigned ints. If m_RefBCO is < m_mvtx_bco_range then we will overflow and delete all hits
  while (m_MvtxRawHitMap.begin()->first < ref_bco_minus_range)
  {
 
    if (Verbosity() > 2)
    {
      std::cout << "ditching mvtx bco 0x" << std::hex << m_MvtxRawHitMap.begin()->first << ", ref: 0x" << m_RefBCO << std::dec << std::endl;
    }
    for (auto iter : m_MvtxInputVector)
    {
      iter->CleanupUsedPackets(m_MvtxRawHitMap.begin()->first);
      iter->clearFeeGTML1BCOMap(m_MvtxRawHitMap.begin()->first);
    }
    for (auto mvtxFeeIdInfo : m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector)
    {
      if (Verbosity() > 1)
      {
        mvtxFeeIdInfo->identify();
      }
      delete mvtxFeeIdInfo;
    }
    
    m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector.clear();
    m_MvtxRawHitMap.begin()->second.MvtxL1TrgBco.clear();
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
  

  unsigned int refbcobitshift = m_RefBCO & 0x3FU;
  h_refbco_mvtx->Fill(refbcobitshift);
  for (auto &[strbbco, mvtxrawhitinfo] : m_MvtxRawHitMap)
  {
    auto diff = (m_RefBCO > strbbco) ? m_RefBCO - strbbco : strbbco - m_RefBCO;
    if (diff > m_mvtx_bco_range)
    {
      continue;
    }
    if (diff > (m_RefBCO + m_mvtx_bco_range))
    {
      break;
    }
    for (auto feeidinfo : mvtxrawhitinfo.MvtxFeeIdInfoVector)
    {
      auto feeId = feeidinfo->get_feeId();

      auto link = MvtxRawDefs::decode_feeid(feeId);
      auto [felix, endpoint] = MvtxRawDefs::get_flx_endpoint(link.layer, link.stave);
      int packetid = felix * 2 + endpoint;
      h_tagStBcoFelix_mvtx[packetid]->Fill(refbcobitshift);
      h_tagStBcoFEE_mvtx->Fill(feeId);
    }
    break;
  }

  std::map<int, std::set<int>> taggedPacketsFEEs;
  for (auto &p : m_MvtxInputVector)
  {
    auto gtml1bcoset_perfee = p->getFeeGTML1BCOMap();
//    int feecounter = 0;
    for (auto &[feeid, gtmbcoset] : gtml1bcoset_perfee)
    {
      auto link = MvtxRawDefs::decode_feeid(feeid);
      auto [felix, endpoint] = MvtxRawDefs::get_flx_endpoint(link.layer, link.stave);
      int packetid = felix * 2 + endpoint;
      for (auto &gtmbco : gtmbcoset)
      {
        auto diff = (m_RefBCO > gtmbco) ? m_RefBCO - gtmbco : gtmbco - m_RefBCO;

        h_bcoGL1LL1diff[packetid]->Fill(diff);

        if (diff < 3)
        {
          taggedPacketsFEEs[packetid].insert(feeid);
          break;
        }
      }
//      feecounter++;
    }

  }
  int allfeestagged = 0;
  for (auto &[pid, feeset] : taggedPacketsFEEs)
  {
    h_tagBcoFelix_mvtx[pid]->Fill(refbcobitshift);
    if (feeset.size() == 12)
    {
      allfeestagged++;
      h_tagBcoFelixAllFees_mvtx[pid]->Fill(refbcobitshift);
    }
    feeset.clear();
  }
  if (allfeestagged == 12)
  {
    h_taggedAllFelixesAllFees_mvtx->Fill(refbcobitshift);
  }
  if (taggedPacketsFEEs.size() == 12)
  {
    h_taggedAllFelixes_mvtx->Fill(refbcobitshift);
  }
  taggedPacketsFEEs.clear();
  
  if (m_mvtx_is_triggered)
  {
    while (select_crossings <= m_MvtxRawHitMap.begin()->first && m_MvtxRawHitMap.begin()->first <= select_crossings + m_mvtx_bco_range) //triggered
    {
      if (Verbosity() > 2)
      {
        std::cout << "Adding 0x" << std::hex << m_MvtxRawHitMap.begin()->first
                  << " ref: 0x" << select_crossings << std::dec << std::endl;
      }
      for (auto mvtxFeeIdInfo : m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector)
      {
        if (Verbosity() > 1)
        {
          mvtxFeeIdInfo->identify();
        }
        mvtxEvtHeader->AddFeeIdInfo(mvtxFeeIdInfo);
        delete mvtxFeeIdInfo;
      }
      m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector.clear();
      mvtxEvtHeader->AddL1Trg(m_MvtxRawHitMap.begin()->second.MvtxL1TrgBco);

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
      m_MvtxRawHitMap.begin()->second.MvtxL1TrgBco.clear();
      m_MvtxRawHitMap.erase(m_MvtxRawHitMap.begin());
      // m_MvtxRawHitMap.empty() need to be checked here since we do not call FillPoolMvtx()
      if (m_MvtxRawHitMap.empty())
      {
        break;
      }
    }
  }
  else
  {
    while (m_MvtxRawHitMap.begin()->first <= select_crossings - m_mvtx_bco_range) //streamed
    {
      if (Verbosity() > 2)
      {
        std::cout << "Adding 0x" << std::hex << m_MvtxRawHitMap.begin()->first
                  << " ref: 0x" << select_crossings << std::dec << std::endl;
      }
      for (auto mvtxFeeIdInfo : m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector)
      {
        if (Verbosity() > 1)
        {
          mvtxFeeIdInfo->identify();
        }
        mvtxEvtHeader->AddFeeIdInfo(mvtxFeeIdInfo);
        delete mvtxFeeIdInfo;
      }
      m_MvtxRawHitMap.begin()->second.MvtxFeeIdInfoVector.clear();
      mvtxEvtHeader->AddL1Trg(m_MvtxRawHitMap.begin()->second.MvtxL1TrgBco);

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
      m_MvtxRawHitMap.begin()->second.MvtxL1TrgBco.clear();
      m_MvtxRawHitMap.erase(m_MvtxRawHitMap.begin());
      // m_MvtxRawHitMap.empty() need to be checked here since we do not call FillPoolMvtx()
      if (m_MvtxRawHitMap.empty())
      {
        break;
      }
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
  if (!container)
  {
    container = findNode::getClass<MicromegasRawHitContainer>(m_topNode, (*(m_MicromegasInputVector.begin()))->getHitContainerName());
    if (!container)
    {
      std::cout << PHWHERE << "No micromegas raw hit container found, exiting." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  // get reference BCO from Micromegas data stream if not set already
  if (m_RefBCO == 0)
  {
    m_RefBCO = m_MicromegasRawHitMap.begin()->first;
  }

  // define bco range
  const uint64_t first_bco = m_RefBCO - m_micromegas_negative_bco;
  const uint64_t last_bco = m_RefBCO + m_micromegas_bco_range - m_micromegas_negative_bco;
  if (Verbosity() > 2)
  {
    std::cout
      << "Fun4AllStreamingInputManager::FillMicromegas - select Micromegas crossings"
      << " from 0x" << std::hex << first_bco
      << " to 0x" << last_bco
      << " for ref BCO " << m_RefBCO
      << std::dec << std::endl;
  }

  // cleanup all data that correspond to too early BCO. Said data is effectively dropped
  while (m_MicromegasRawHitMap.begin()->first < first_bco)
  {
    if (Verbosity() > 2)
    {
      std::cout
        << "Micromegas BCO: 0x" << std::hex << m_MicromegasRawHitMap.begin()->first
        << " smaller than GL1 BCO: 0x" << first_bco
        << ", ditching this bco" << std::dec << std::endl;
    }

    for (const auto& poolinput : m_MicromegasInputVector)
    { poolinput->CleanupUsedPackets(m_MicromegasRawHitMap.begin()->first, true); }

    // remove
    m_MicromegasRawHitMap.erase(m_MicromegasRawHitMap.begin());

    // fill pools again
    iret = FillMicromegasPool();
    if (iret)
    {
      return iret;
    }
  }

  // fill all BCO statistics
  for (const auto &iter : m_MicromegasInputVector)
  {
    iter->FillBcoQA(m_RefBCO);
  }

  // store hits relevant for this trigger and cleanup
  for( auto iter = m_MicromegasRawHitMap.begin(); iter != m_MicromegasRawHitMap.end() && iter->first <= last_bco;  iter = m_MicromegasRawHitMap.erase(iter))
  {
    for (const auto &hititer : iter->second.MicromegasRawHitVector)
    {
      container->AddHit(hititer);
    }

    for (const auto &poolinput : m_MicromegasInputVector)
    {
      poolinput->CleanupUsedPackets(iter->first);
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
  if (!tpccont)
  {
    /// if we set the node name and are running over single prdfs, thre is only one prdf in the vector
    tpccont = findNode::getClass<TpcRawHitContainer>(m_topNode, (*(m_TpcInputVector.begin()))->getHitContainerName());
    if (!tpccont)
    {
      std::cout << PHWHERE << "No tpc raw hit container found, exiting." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }
  //  std::cout << "before filling m_TpcRawHitMap size: " <<  m_TpcRawHitMap.size() << std::endl;
  uint64_t select_crossings = m_tpc_bco_range;
  if (m_RefBCO == 0)
  {
    m_RefBCO = m_TpcRawHitMap.begin()->first;
  }
  select_crossings += m_RefBCO;
  if (Verbosity() > 2)
  {

    std::cout << "select TPC crossings"
              << " from 0x" << std::hex << m_RefBCO - m_tpc_negative_bco
              << " to 0x" << select_crossings - m_tpc_negative_bco
              << std::dec << std::endl;

  }
  // m_TpcRawHitMap.empty() does not need to be checked here, FillTpcPool returns non zero
  // if this map is empty which is handled above

  if (m_TpcRawHitMap.size()>0)
  {
    while (m_TpcRawHitMap.begin()->first < m_RefBCO - m_tpc_negative_bco)
    {
      for (auto iter : m_TpcInputVector)
      {
        iter->CleanupUsedPackets(m_TpcRawHitMap.begin()->first);
        iter->clearPacketBClkStackMap(m_TpcRawHitMap.begin()->first);

    }
      m_TpcRawHitMap.begin()->second.TpcRawHitVector.clear();
      m_TpcRawHitMap.erase(m_TpcRawHitMap.begin());
      iret = FillTpcPool();
      if (iret)
      {
        return iret;
      }
    }
  }
  unsigned int refbcobitshift = m_RefBCO & 0x3FU;
  h_refbco_tpc->Fill(refbcobitshift);
  bool allpackets = true;
  for (auto &p : m_TpcInputVector)
  {
    auto bcl_stack = p->BclkStackMap();
    int packetnum = 0;
    int histo_to_fill = (bcl_stack.begin()->first - 4000) / 10;

    for (auto &[packetid, bclset] : bcl_stack)
    {
      bool thispacket = false;
      for (auto &bcl : bclset)
      {
        auto diff = (m_RefBCO > bcl) ? m_RefBCO - bcl : bcl - m_RefBCO;
        if (diff < 256)
        {
          thispacket = true;
          h_gl1tagged_tpc[histo_to_fill][packetnum]->Fill(refbcobitshift);
        }
      }
      if (thispacket == false)
      {
        allpackets = false;
      }
    
      packetnum++;
    }
  }
  if (allpackets)
  {
    h_taggedAll_tpc->Fill(refbcobitshift);
  }
  // again m_TpcRawHitMap.empty() is handled by return of FillTpcPool()
  if (m_TpcRawHitMap.size()>0)
  {
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
          // we just want to erase anything that is well away from the current GL1
  
    }
      m_TpcRawHitMap.begin()->second.TpcRawHitVector.clear();
      m_TpcRawHitMap.erase(m_TpcRawHitMap.begin());
      if (m_TpcRawHitMap.empty())
      {
        break;
      }
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
  uint64_t ref_bco_minus_range = 0;
  if (m_RefBCO > m_intt_negative_bco)
  {
    ref_bco_minus_range = m_RefBCO - m_intt_negative_bco;
  }
  for (auto iter : m_InttInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillInttPool - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(ref_bco_minus_range);
    // iter->FillPool();
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
    std::cout << "InttRawHitMap is empty - we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillTpcPool()
{
  uint64_t ref_bco_minus_range = 0;
  if(m_RefBCO > m_tpc_negative_bco)
  {
    ref_bco_minus_range = m_RefBCO - m_tpc_negative_bco;
  }

  for (auto iter : m_TpcInputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingInputManager::FillTpcPool - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(ref_bco_minus_range);
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
  // if (m_TpcRawHitMap.empty())
  // {
  //   std::cout << "TpcRawHitMap is empty - we are done" << std::endl;
    // return -1;
  // }
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
    std::cout << "MicromegasRawHitMap is empty - we are done" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllStreamingInputManager::FillMvtxPool()
{
  uint64_t ref_bco_minus_range = m_RefBCO < m_mvtx_bco_range ? 0 : m_RefBCO - m_mvtx_bco_range;
  for (auto iter : m_MvtxInputVector)
  {
    if (Verbosity() > 3)
    {
      std::cout << "Fun4AllStreamingInputManager::FillMvtxPool - fill pool for " << iter->Name() << std::endl;
    }
    iter->FillPool(ref_bco_minus_range);
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
    
    std::cout << "MvtxRawHitMap is empty - we are done" << std::endl;
    return -1;
  }
  return 0;
}
void Fun4AllStreamingInputManager::createQAHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH1I("h_TpcPoolQA_RefGL1BCO", "TPC ref BCO", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_InttPoolQA_RefGL1BCO", "INTT ref BCO", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_MvtxPoolQA_RefGL1BCO", "MVTX ref BCO", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }

  {
    auto h = new TH1I("h_InttPoolQA_TagBCOAllServers", "INTT trigger tagged BCO all servers", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_MvtxPoolQA_TagBCOAllFelixs", "MVTX trigger tagged BCO all felixs", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_MvtxPoolQA_TagBCOAllFelixsAllFees", "MVTX trigger tagged BCO all felixes and fees", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    hm->registerHisto(h);
    h_taggedAllFelixesAllFees_mvtx = h;
  }
  {
    auto h = new TH1I("h_TpcPoolQA_TagBCOAllPackets", "TPC trigger tagged BCO all servers", 1000,
                      0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h);
  }
  h_tagStBcoFEE_mvtx = new TH1I("h_MvtxPoolQA_TagStBcoFEEs", "", 10000, 0, 10000);
  hm->registerHisto(h_tagStBcoFEE_mvtx);

  // intt has 8 prdfs, one per felix
  for (int i = 0; i < 8; i++)
  {
    auto h = new TH1I((boost::format("h_InttPoolQA_TagBCO_server%i") % i).str().c_str(), "INTT trigger tagged BCO", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("GL1 BCO");
    h->SetTitle((boost::format("EBDC %i") % i).str().c_str());
    hm->registerHisto(h);

    auto h_all = new TH1I((boost::format("h_InttPoolQA_TagBCOAllFees_Server%i") % i).str().c_str(), "INTT trigger tagged BCO all servers", 1000, 0, 1000);
    h_all->GetXaxis()->SetTitle("GL1 BCO");
    h_all->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h_all);
    for (int j = 0; j < 14; j++)
    {
      auto h2 = new TH1I((boost::format("h_InttPoolQA_TagBCO_server%i_fee%i") % i % j).str().c_str(), "INTT trigger tagged BCO per FEE", 1000, 0, 1000);
      h2->GetXaxis()->SetTitle("GL1 BCO");
      h2->SetTitle((boost::format("EBDC %i FEE %i") % i % j).str().c_str());
      hm->registerHisto(h2);
    }
  }
  for (int i = 0; i < 12; i++)
  {
    {
      auto h = new TH1I((boost::format("h_MvtxPoolQA_TagBCO_felix%i") % i).str().c_str(), "MVTX trigger tagged BCO", 1000, 0, 1000);
      h->GetXaxis()->SetTitle("GL1 BCO");
      h->SetTitle((boost::format("Felix %i") % i).str().c_str());
      hm->registerHisto(h);
    }
    {
      auto h = new TH1I((boost::format("h_MvtxPoolQA_TagStBco_felix%i") % i).str().c_str(), "", 1000, 0, 1000);
      hm->registerHisto(h);
      h_tagStBcoFelix_mvtx[i] = h;
    }
    auto h_all = new TH1I((boost::format("h_MvtxPoolQA_TagBCOAllFees_Felix%i") % i).str().c_str(), "MVTX trigger tagged BCO all Fees", 1000, 0, 1000);
    h_all->GetXaxis()->SetTitle("GL1 BCO");
    h_all->SetTitle("GL1 Reference BCO");
    hm->registerHisto(h_all);
  }
  for (int i = 0; i < 24; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      {
        auto h = new TH1I((boost::format("h_TpcPoolQA_TagBCO_ebdc%i_packet%i") % i % j).str().c_str(), "TPC trigger tagged BCO", 1000, 0, 1000);
        h->GetXaxis()->SetTitle("GL1 BCO");
        h->SetTitle((boost::format("Packet %i and packet %i") % i % j).str().c_str());
        hm->registerHisto(h);
      }
    }
  }

  for (int i = 0; i < 12; i++)
  {
    h_bcoGL1LL1diff[i] = new TH1I((boost::format("h_MvtxPoolQA_GL1LL1BCODiff_packet%i") % i).str().c_str(), "MVTX BCO diff;|GL1 BCO - LL1 BCO|", 5000, 0, 5000);
    hm->registerHisto(h_bcoGL1LL1diff[i]);
    h_bcoLL1Strobediff[i] = new TH1I((boost::format("h_MvtxPoolQA_LL1StrobeBCODiff_packet%i") % i).str().c_str(), "MVTX BCO diff; |LL1 BCO - Strobe BCO|", 100000, 0, 100000);
    hm->registerHisto(h_bcoLL1Strobediff[i]);
  }
  // Get the global pointers
  h_refbco_intt = dynamic_cast<TH1 *>(hm->getHisto("h_InttPoolQA_RefGL1BCO"));
  h_taggedAll_intt = dynamic_cast<TH1 *>(hm->getHisto("h_InttPoolQA_TagBCOAllServers"));
  h_taggedAllFee_intt = new TH1I("h_InttPoolQA_TagBCOAllServersAllFees", "INTT trigger tagged BCO all servers and fees", 1000, 0, 1000);
  hm->registerHisto(h_taggedAllFee_intt);
  for (int i = 0; i < 8; i++)
  {
    h_gl1tagged_intt[i] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_InttPoolQA_TagBCO_server%i") % i).str().c_str()));
    for (int j = 0; j < 14; j++)
    {
      h_gl1taggedfee_intt[i][j] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_InttPoolQA_TagBCO_server%i_fee%i") % i % j).str().c_str()));
    }
    h_taggedAllFees_intt[i] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_InttPoolQA_TagBCOAllFees_Server%i") % i).str().c_str()));
  }

  h_refbco_mvtx = dynamic_cast<TH1 *>(hm->getHisto("h_MvtxPoolQA_RefGL1BCO"));
  h_taggedAllFelixes_mvtx = dynamic_cast<TH1 *>(hm->getHisto("h_MvtxPoolQA_TagBCOAllFelixs"));
  for (int i = 0; i < 12; i++)
  {
    h_tagBcoFelix_mvtx[i] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_MvtxPoolQA_TagBCO_felix%i") % i).str().c_str()));
    h_tagBcoFelixAllFees_mvtx[i] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_MvtxPoolQA_TagBCOAllFees_Felix%i") % i).str().c_str()));
  }

  for (int i = 0; i < 24; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      h_gl1tagged_tpc[i][j] = dynamic_cast<TH1 *>(hm->getHisto((boost::format("h_TpcPoolQA_TagBCO_ebdc%i_packet%i") % i % j).str()));
    }
  }

  h_refbco_tpc = dynamic_cast<TH1 *>(hm->getHisto("h_TpcPoolQA_RefGL1BCO"));
  h_taggedAll_tpc = dynamic_cast<TH1 *>(hm->getHisto("h_TpcPoolQA_TagBCOAllPackets"));
}

#include "SingleTpcTimeFrameInput.h"
#include "TpcTimeFrameBuilder.h"
#include "TpcTimeFrameBuilderRun3.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <qautils/QAHistManagerDef.h>

#include <ffarawobjects/TpcRawHitContainerv3.h>
#include <ffarawobjects/TpcRawHitv3.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/InputFileHandlerReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/oncsSubConstants.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <memory>
#include <set>

SingleTpcTimeFrameInput::SingleTpcTimeFrameInput(const std::string &name)
  : SingleStreamingInput(name)
  , plist(new Packet *[NTPCPACKETS])
{
  SubsystemEnum(InputManagerType::TPC);
  m_rawHitContainerName = "TPCRAWHIT";

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_FillPoolTimer = new PHTimer("SingleTpcTimeFrameInput_" + name + "_FillPool");
  m_getNextEventTimer = new PHTimer("SingleTpcTimeFrameInput_" + name + "_getNextEvent");
  m_ProcessPacketTimer = new PHTimer("SingleTpcTimeFrameInput_" + name + "_ProcessPacket");
  m_getTimeFrameTimer = new PHTimer("SingleTpcTimeFrameInput_" + name + "_getTimeFrame");

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  m_hNorm = new TH1D(TString("SingleTpcTimeFrameInput_" + name) + "_Normalization",  //
                     TString("SingleTpcTimeFrameInput_" + name) + " Normalization;Items;Count",
                     20, .5, 20.5);
  int i = 1;
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FillPoolCount");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FillPoolTime[ms]");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "getNextEventCount");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "getNextEventTime[ms]");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "ProcessPacketCount");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "ProcessPacketTime[ms]");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "getTimeFrameCount");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "getTimeFrameTime[ms]");
  assert(i <= 20);
  m_hNorm->GetXaxis()->LabelsOption("v");
  hm->registerHisto(m_hNorm);
}

SingleTpcTimeFrameInput::~SingleTpcTimeFrameInput()
{
  delete[] plist;
  delete m_FillPoolTimer;
  delete m_getNextEventTimer;
  delete m_ProcessPacketTimer;
  delete m_getTimeFrameTimer;
  for (auto iter : m_TpcTimeFrameBuilderMap)
  {
    delete iter.second;
  }
}

SingleTpcTimeFrameInput::TimeTracker::TimeTracker(PHTimer *timer, const std::string &name, TH1 *hout)
  : m_timer(timer)
  , m_name(name)
  , m_hNorm(hout)
{
  assert(m_timer);
  assert(m_hNorm);

  m_hNorm->Fill((m_name + "Count").c_str(), 1);

  m_timer->restart();
}

void SingleTpcTimeFrameInput::TimeTracker::stop()
{
  assert(m_timer);
  assert(m_hNorm);

  m_timer->stop();
  m_hNorm->Fill((m_name + "Time[ms]").c_str(), m_timer->elapsed());

  stopped = true;
}

SingleTpcTimeFrameInput::TimeTracker::~TimeTracker()
{
  if (!stopped)
  {
    stop();
  }
}

void SingleTpcTimeFrameInput::FillPool(const uint64_t targetBCO)
{
  m_FillPoolStatus = Fun4AllReturnCodes::EVENT_OK;
  {
    static bool first = true;
    if (first)
    {
      first = false;

      if (!m_SelectedPacketIDs.empty())
      {
        std::cout << "SingleTpcTimeFrameInput::" << Name() << " : note, only processing packets with ID: ";
        for (const auto &id : m_SelectedPacketIDs)
        {
          std::cout << id << " ";
        }
        std::cout << std::endl;
      }
    }
  }

  TimeTracker fillPoolTimer(m_FillPoolTimer, "FillPool", m_hNorm);

  if ((Verbosity() >= 1 && targetBCO % 941 == 10) || Verbosity() >= 2)
  {
    std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
              << " Entry with targetBCO = 0x" << std::hex << targetBCO
              << "(" << std::dec << targetBCO << ")" << std::endl;

    m_FillPoolTimer->print_stat();
    m_getNextEventTimer->print_stat();
    m_ProcessPacketTimer->print_stat();
    m_getTimeFrameTimer->print_stat();
  }

  if (AllDone())  // no more files and all events read
  {
    if (Verbosity() > 1)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " AllDone for targetBCO " << targetBCO << std::endl;
    }

    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (Verbosity() > 1)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " GetEventiterator == null for targetBCO " << targetBCO << std::endl;
    }

    if (OpenNextFile() == InputFileHandlerReturnCodes::FAILURE)
    {
      if (Verbosity() > 1)
      {
        std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                  << " OpenNextFile() quit for targetBCO " << targetBCO << std::endl;
      }
      AllDone(1);
      return;
    }
  }
  //  std::set<uint64_t> saved_beamclocks;
  while (true)
  {
    // clean up cache to avoid memory over usage when trigger jumped by a long time
    CleanupUsedPackets(targetBCO - kUsedPacketsCachingLimit);

    if (m_TpcTimeFrameBuilderMap.empty())
    {
      if (Verbosity() > 1)
      {
        std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                  << " m_TpcTimeFrameBuilderMap empty for targetBCO 0x"
                  << std::hex << targetBCO << std::dec << ". Start processing next event... " << std::endl;
      }
    }
    else
    {
      bool require_more_data = false;
      for (auto &map_builder : m_TpcTimeFrameBuilderMap)
      {
        require_more_data |= map_builder.second->isMoreDataRequired(targetBCO);
      }
      if (!require_more_data)
      {
        if (Verbosity() > 1)
        {
          std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                    << " satisified require_more_data for targetBCO 0x"
                    << std::hex << targetBCO << std::dec << std::endl;
        }

        break;
      }
    }
    if (Verbosity() > 1)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " require_more_data for targetBCO 0x"
                << std::hex << targetBCO << std::dec << std::endl;
    }

    TimeTracker getNextEventTimer(m_getNextEventTimer, "getNextEvent", m_hNorm);
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (OpenNextFile() == InputFileHandlerReturnCodes::FAILURE)
      {
        AllDone(1);
        return;
      }
      evt.reset(GetEventiterator()->getNextEvent());
      RunNumber(0);
    }

    if (RunNumber() == 0)
    {
      RunNumber(evt->getRunNumber());

      if (Verbosity() >= 1)
      {
        std::cout << __PRETTY_FUNCTION__ << ": Fetching new file " << FileName()
                  << " for run " << evt->getRunNumber() << " with next Event # " << evt->getEvtSequence() << std::endl;
      }
    }

    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }
    // int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, NTPCPACKETS);

    if (npackets >= NTPCPACKETS)
    {
      std::cout << PHWHERE << " Packets array size " << NTPCPACKETS
                << " too small for " << Name()
                << ", increase NTPCPACKETS and rebuild" << std::endl;
      exit(1);
    }
    getNextEventTimer.stop();

    TimeTracker ProcessPacketTimer(m_ProcessPacketTimer, "ProcessPacket", m_hNorm);
    for (int i = 0; i < npackets; i++)
    {
      // keep pointer to local packet
      auto &packet = plist[i];
      assert(packet);

      auto cleanup_remaining_packets = [&](const int first_index)
      {
        for (int j = first_index; j < npackets; ++j)
        {
          if (plist[j])
          {
            delete plist[j];
            plist[j] = nullptr;
          }
        }
      };

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (!m_SelectedPacketIDs.empty() && !m_SelectedPacketIDs.contains(packet_id))
      {
        if (Verbosity() > 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": Skipping packet id: " << packet_id << std::endl;
        }

        delete packet;
        packet = nullptr;
        continue;
      }

      const int hit_format = packet->getHitFormat();
      const auto builder_hit_format_iter = m_TpcTimeFrameBuilderHitFormatMap.find(packet_id);
      if (builder_hit_format_iter != m_TpcTimeFrameBuilderHitFormatMap.end() && builder_hit_format_iter->second != hit_format)
      {
        std::cout << __PRETTY_FUNCTION__ << ": Error : packet id " << packet_id
                  << " changed TPC hit format from " << builder_hit_format_iter->second
                  << " to " << hit_format << ". Aborting run." << std::endl;
        packet->identify();
        m_FillPoolStatus = Fun4AllReturnCodes::ABORTRUN;
        cleanup_remaining_packets(i);
        return;
      }

      if (!m_TpcTimeFrameBuilderMap.contains(packet_id))
      {
        TpcTimeFrameBuilderBase *builder = nullptr;
        if (hit_format == IDTPCFEEV4)
        {
          if (Verbosity() >= 1)
          {
            std::cout << __PRETTY_FUNCTION__ << ": Creating TpcTimeFrameBuilder for packet id: " << packet_id
                      << " hit format " << hit_format << std::endl;
          }
          builder = new TpcTimeFrameBuilder(packet_id);
        }
        else if (hit_format == IDTPCFEEV5 || hit_format == IDTPCFEEV6)
        {
          if (Verbosity() >= 1)
          {
            std::cout << __PRETTY_FUNCTION__ << ": Creating TpcTimeFrameBuilderRun3 for packet id: " << packet_id
                      << " hit format " << hit_format << std::endl;
          }
          builder = new TpcTimeFrameBuilderRun3(packet_id);
        }
        else
        {
          std::cout << __PRETTY_FUNCTION__ << ": Error : unsupported TPC hit format " << hit_format
                    << " for packet id " << packet_id << ". Aborting run." << std::endl;
          packet->identify();
          m_FillPoolStatus = Fun4AllReturnCodes::ABORTRUN;
          cleanup_remaining_packets(i);
          return;
        }

        m_TpcTimeFrameBuilderMap[packet_id] = builder;
        m_TpcTimeFrameBuilderHitFormatMap[packet_id] = hit_format;
        m_TpcTimeFrameBuilderMap[packet_id]->setVerbosity(Verbosity());
        m_TpcTimeFrameBuilderMap[packet_id]->fillBadFeeMap();
        if (!m_digitalCurrentDebugTTreeName.empty())
        {
          m_TpcTimeFrameBuilderMap[packet_id]->SaveDigitalCurrentDebugTTree(m_digitalCurrentDebugTTreeName);
        }
        if (!m_bxCounterSyncCDBTTreeName.empty())
        {
          m_TpcTimeFrameBuilderMap[packet_id]->SaveBXCounterSyncCDBTTree(m_bxCounterSyncCDBTTreeName);
        }
      }

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      assert(m_TpcTimeFrameBuilderMap[packet_id]);
      const int process_packet_status = m_TpcTimeFrameBuilderMap[packet_id]->ProcessPacket(packet);
      if (process_packet_status < 0)
      {
        std::cout << __PRETTY_FUNCTION__ << ": Error : TPC packet builder returned " << process_packet_status
                  << " for packet id " << packet_id << ". Aborting run." << std::endl;
        m_FillPoolStatus = process_packet_status;
        cleanup_remaining_packets(i);
        return;
      }
      // require_more_data = require_more_data or m_TpcTimeFrameBuilderMap[packet_id]->isMoreDataRequired(targetBCO);

      delete packet;
      packet = nullptr;
    }  //     for (int i = 0; i < npackets; i++)
    ProcessPacketTimer.stop();

  }  // while (require_more_data)

  // output the time frame
  for (auto &map_builder : m_TpcTimeFrameBuilderMap)
  {
    assert(!map_builder.second->isMoreDataRequired(targetBCO));
    auto &timeframe = map_builder.second->getTimeFrame(targetBCO);

    for (auto *newhit : timeframe)
    {
      StreamingInputManager()->AddTpcRawHit(targetBCO, newhit);
    }
  }

  TimeTracker getTimeFrameTimer(m_getTimeFrameTimer, "getTimeFrame", m_hNorm);
}

void SingleTpcTimeFrameInput::Print(const std::string & /*what*/) const
{
}

void SingleTpcTimeFrameInput::CleanupUsedPackets(const uint64_t bclk)
{
  if (Verbosity() > 2)
  {
    std::cout << "cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }

  for (auto &iter : m_TpcTimeFrameBuilderMap)
  {
    iter.second->CleanupUsedPackets(bclk);
  }
}

void SingleTpcTimeFrameInput::ClearCurrentEvent()
{
  std::cout << "SingleTpcTimeFrameInput::ClearCurrentEvent() - deprecated " << std::endl;
  exit(1);

  // // called interactively, to get rid of the current event
  // uint64_t currentbclk = *m_BclkStack.begin();
  // //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  // CleanupUsedPackets(currentbclk);
  // // m_BclkStack.erase(currentbclk);
  // // m_BeamClockFEE.erase(currentbclk);
  return;
}

void SingleTpcTimeFrameInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "TPC"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("TPC");
    dstNode->addNode(detNode);
  }
  TpcRawHitContainer *tpchitcont = findNode::getClass<TpcRawHitContainer>(detNode, m_rawHitContainerName);
  if (!tpchitcont)
  {
    tpchitcont = new TpcRawHitContainerv3();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(tpchitcont, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

void SingleTpcTimeFrameInput::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetTpcBcoRange(m_BcoRange);
    StreamingInputManager()->SetTpcNegativeBco(m_NegativeBco);
  }
  return;
}

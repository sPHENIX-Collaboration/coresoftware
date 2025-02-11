#include "SingleTpcTimeFrameInput.h"
#include "TpcTimeFrameBuilder.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/TpcRawHitContainerv3.h>
#include <ffarawobjects/TpcRawHitv3.h>

#include <frog/FROG.h>
#include <phool/PHTimer.h>  // for PHTimer

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <qautils/QAHistManagerDef.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <memory>
#include <set>

SingleTpcTimeFrameInput::SingleTpcTimeFrameInput(const std::string &name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::TPC);
  m_rawHitContainerName = "TPCRAWHIT";
  plist = new Packet *[NTPCPACKETS];

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
  , m_name(name.c_str())
  , m_hNorm(hout)
{
  assert(m_timer);
  assert(m_hNorm);

  m_hNorm->Fill(m_name + "Count", 1);

  m_timer->restart();
}

void SingleTpcTimeFrameInput::TimeTracker::stop()
{
  assert(m_timer);
  assert(m_hNorm);

  m_timer->stop();
  m_hNorm->Fill(m_name + "Time[ms]", m_timer->elapsed());

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
  {
    static bool first = true;
    if (first)
    {
      first = false;

      if (m_SelectedPacketIDs.size())
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

  if ((Verbosity() >= 1 and targetBCO % 941 == 10) or Verbosity() >= 2)
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

    if (!OpenNextFile())
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
      if (not require_more_data)
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
      if (!OpenNextFile())
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
        std::cout << __PRETTY_FUNCTION__ << ": Fetching new file "<<FileName()
        <<" for run "<<evt->getRunNumber()<<" with next Event # " << evt->getEvtSequence() << std::endl;
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

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (m_SelectedPacketIDs.size() > 0 and m_SelectedPacketIDs.find(packet_id) == m_SelectedPacketIDs.end())
      {
        if (Verbosity() > 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": Skipping packet id: " << packet_id << std::endl;
        }

        delete packet;
        packet = nullptr;
        continue;
      }

      if (m_TpcTimeFrameBuilderMap.find(packet_id) == m_TpcTimeFrameBuilderMap.end())
      {
        if (Verbosity() >= 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": Creating TpcTimeFrameBuilder for packet id: " << packet_id << std::endl;
        }

        m_TpcTimeFrameBuilderMap[packet_id] = new TpcTimeFrameBuilder(packet_id);
        m_TpcTimeFrameBuilderMap[packet_id]->setVerbosity(Verbosity());
      }

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      assert(m_TpcTimeFrameBuilderMap[packet_id]);
      m_TpcTimeFrameBuilderMap[packet_id]->ProcessPacket(packet);
      // require_more_data = require_more_data or m_TpcTimeFrameBuilderMap[packet_id]->isMoreDataRequired(targetBCO);

      delete packet;
      packet = nullptr;
    }  //     for (int i = 0; i < npackets; i++)
    ProcessPacketTimer.stop();

  }  // while (require_more_data)

  // output the time frame
  for (auto &map_builder : m_TpcTimeFrameBuilderMap)
  {
    assert(not map_builder.second->isMoreDataRequired(targetBCO));
    auto &timeframe = map_builder.second->getTimeFrame(targetBCO);

    for (auto newhit : timeframe)
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

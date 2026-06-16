#include "TpcTimeFrameBuilderRun3.h"

#include <Event/oncsSubConstants.h>
#include <Event/packet.h>

#include <qautils/QAHistManagerDef.h>

#include <ffarawobjects/TpcRawHitv2.h>
#include <ffarawobjects/TpcRawHitv3.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <phool/PHTimer.h>  // for PHTimer

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TTree.h>
#include <TVector3.h>

#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <tuple>  // For std::tie

TpcTimeFrameBuilderRun3::TpcTimeFrameBuilderRun3(const int packet_id)
  : m_packet_id(packet_id)
  , m_HistoPrefix("TpcTimeFrameBuilderRun3_Packet" + std::to_string(packet_id))
  , m_bxCounterSyncCDBTTreeName(m_HistoPrefix + "_BXCounterSyncCDBTTree.root")
{
  for (int fee = 0; fee < MAX_FEECOUNT; ++fee)
  {
    m_bcoMatchingInformation_vec.emplace_back(
        std::string("BcoMatchingInformation_Packet") + std::to_string(packet_id) + "_FEE" + std::to_string(fee));
  }

  m_feeData.resize(MAX_FEECOUNT);
  m_timeHitMap.resize(MAX_FEECOUNT);

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_packetTimer = new PHTimer("TpcTimeFrameBuilderRun3_Packet" + std::to_string(packet_id));

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  m_hNorm = new TH1D(TString(m_HistoPrefix.c_str()) + "_Normalization",  //
                     TString(m_HistoPrefix.c_str()) + " Normalization;Items;Count",
                      kRun3NormalizationBinCount, .5, kRun3NormalizationBinCount + .5);
  int i = 1;
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Packet");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Lv1-Taggers");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "EnDat-Taggers");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "ChannelPackets");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Waveforms");

  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_GTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_FEE");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_FEE_INVALID");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_INVALID");

  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_GTM_HEARTBEAT");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "DMA_WORD_GTM_DC_STOP_SEND");

  m_hNorm->GetXaxis()->SetBinLabel(i++, "TimeFrameSizeLimitError");

  m_hNorm->GetXaxis()->SetBinLabel(i++, "GTM_TimeFrame_Matched");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "GTM_TimeFrame_Unmatched");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "GTM_TimeFrame_Matched_Hit_Sum");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "GTM_TimeFrame_Dropped_Hit_Sum");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Run3_TimeFrame_Exact_Matched");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Run3_TimeFrame_FuzzyFallback");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Run3_TimeFrame_FuzzyFallback_Hit_Sum");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "Run3_TimeFrame_MatchFailed");

  m_hNormTruncatedWaveformRecoveryFeeFirstBin = i;
  for (uint16_t fee = 0; fee < MAX_FEECOUNT; ++fee)
  {
    m_hNorm->GetXaxis()->SetBinLabel(i++, ("Run3_TruncatedWaveformRecover_FEE" + std::to_string(fee)).c_str());
  }

  assert(i <= kRun3NormalizationBinCount + 1);
  m_hNorm->GetXaxis()->LabelsOption("v");
  hm->registerHisto(m_hNorm);

  h_PacketLength = new TH1I(TString(m_HistoPrefix.c_str()) + "_PacketLength",  //
                            TString(m_HistoPrefix.c_str()) + " PacketLength;PacketLength [32bit Words];Count", 1000, .5, 5e6);
  hm->registerHisto(h_PacketLength);

  h_PacketLength_Residual = new TH1I(TString(m_HistoPrefix.c_str()) + "_PacketLength_Residual",  //
                                     TString(m_HistoPrefix.c_str()) +
                                         " PacketLength that does not fit into DMA transfer;PacketLength [16bit Words];Count",
                                     16, -.5, 15.5);
  hm->registerHisto(h_PacketLength_Residual);

  h_PacketLength_Padding = new TH1I(TString(m_HistoPrefix.c_str()) + "_PacketLength_Padding",  //
                                    TString(m_HistoPrefix.c_str()) +
                                        " padding within PacketLength;PacketLength [32bit Words];Count",
                                    16, -.5, 15.5);
  hm->registerHisto(h_PacketLength_Padding);

  m_hFEEDataStream = new TH2I(TString(m_HistoPrefix.c_str()) + "_FEE_DataStream_WordCount",  //
                              TString(m_HistoPrefix.c_str()) +
                                  " FEE Data Stream Word Count;FEE ID;Type;Count",
                              MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5, 25, .5, 25.5);
  i = 1;
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordValid");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordSkipped");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordDigitalCurrentKeyWord");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "InvalidLength");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "RawHit");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitFormatErrorOverLength");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitFormatErrorMismatchedLength");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitCRCError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "DigitalCurrent");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "DigitalCurrentFormatErrorMismatchedLength");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "DigitalCurrentCRCError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "ParityError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitUnusedBeforeCleanup");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeat");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncUnavailable");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncOK");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncUnavailable");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncOK");
  assert(i <= 25);
  hm->registerHisto(m_hFEEDataStream);

  m_hFEEChannelPacketCount = new TH1I(TString(m_HistoPrefix.c_str()) + "_FEEChannelPacketCount",  //
                                      TString(m_HistoPrefix.c_str()) +
                                          " Count of waveform packet per channel;FEE*256 + Channel;Count",
                                      MAX_FEECOUNT * MAX_CHANNELS, -.5, MAX_FEECOUNT * MAX_CHANNELS - .5);
  hm->registerHisto(m_hFEEChannelPacketCount);

  m_hFEESAMPAADC = new TH2I(TString(m_HistoPrefix.c_str()) + "_FEE_SAMPA_ADC",  //
                            TString(m_HistoPrefix.c_str()) +
                                " ADC distribution in 2D;ADC Time Bin [0...1023];FEE*8+SAMPA;Sum ADC",
                            MAX_PACKET_LENGTH, -.5, MAX_PACKET_LENGTH - .5,
                            MAX_FEECOUNT * MAX_SAMPA, -.5, MAX_FEECOUNT * MAX_SAMPA - .5);
  hm->registerHisto(m_hFEESAMPAADC);

  m_hFEESAMPAHeartBeatSync = new TH1I(TString(m_HistoPrefix.c_str()) + "_FEE_SAMPA_HEARTBEAT_SYNC",  //
                                      TString(m_HistoPrefix.c_str()) +
                                          " FEE/SAMPA Sync Heartbeat Count;FEE*8+SAMPA;Sync Heartbeat Count",
                                      MAX_FEECOUNT * MAX_SAMPA, -.5, MAX_FEECOUNT * MAX_SAMPA - .5);
  hm->registerHisto(m_hFEESAMPAHeartBeatSync);

  h_GTMClockDiff_Matched = new TH1I(TString(m_HistoPrefix.c_str()) + "_GTMClockDiff_Matched",  //
                                    TString(m_HistoPrefix.c_str()) +
                                        " GTM BCO Diff for Matched Time Frame;Trigger BCO Diff [BCO];Count",
                                    1024, -512 - .5, 512 - .5);
  hm->registerHisto(h_GTMClockDiff_Matched);
  h_GTMClockDiff_Unmatched = new TH1I(TString(m_HistoPrefix.c_str()) + "_GTMClockDiff_Unmatched",  //
                                      TString(m_HistoPrefix.c_str()) +
                                          " GTM BCO Diff for Unmatched Time Frame;Trigger BCO Diff [BCO];Count",
                                      1024, -512 - .5, 512 - .5);
  hm->registerHisto(h_GTMClockDiff_Unmatched);
  h_GTMClockDiff_Dropped = new TH1I(TString(m_HistoPrefix.c_str()) + "_GTMClockDiff_Dropped",  //
                                    TString(m_HistoPrefix.c_str()) +
                                        " GTM BCO Diff for Dropped Time Frame;Trigger BCO Diff [BCO];Count",
                                    16384, -16384 - .5, 0 - .5);
  hm->registerHisto(h_GTMClockDiff_Dropped);
  h_TimeFrame_Matched_Size = new TH1I(TString(m_HistoPrefix.c_str()) + "_TimeFrame_Matched_Size",  //
                                      TString(m_HistoPrefix.c_str()) +
                                          " Time frame size for Matched Time Frame ;Size [TPC raw hits];Count",
                                      3328, -.5, 3328 - .5);
  hm->registerHisto(h_TimeFrame_Matched_Size);

  h_Run3_FEE_GTMMatching_ClockDiff = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3_FEE_GTMMatching_ClockDiff",  //
                                                   TString(m_HistoPrefix.c_str()) +
                                                       " Run3 FEE GTM matching clock diff by FEE;Clock Difference [FEE Clock Cycle];FEE;Matched hits",
                                                   2048, -1024 - .5, 1024 - .5,
                                                   MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3_FEE_GTMMatching_ClockDiff);

  h_Run3TimeFrameExactHit_FEE = new TH1I(TString(m_HistoPrefix.c_str()) + "_Run3TimeFrameExactHit_FEE",  //
                                          TString(m_HistoPrefix.c_str()) +
                                              " Run3 exact matched hit sum by FEE;FEE;Exact matched hits",
                                          MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3TimeFrameExactHit_FEE);

  h_Run3TimeFrameFuzzyHit_FEE = new TH1I(TString(m_HistoPrefix.c_str()) + "_Run3TimeFrameFuzzyHit_FEE",  //
                                          TString(m_HistoPrefix.c_str()) +
                                              " Run3 fuzzy fallback matched hit sum by FEE;FEE;Fuzzy fallback matched hits",
                                          MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3TimeFrameFuzzyHit_FEE);

  static constexpr uint32_t kRun3TruncatedWaveformRecoveryWindowPlotingRange = 1200U;
  h_Run3Waveform_GL1Spacing = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3Waveform_GL1Spacing",  //
                                          TString(m_HistoPrefix.c_str()) +
                                              " Run3 matched waveform ADC sum before truncated waveform recovery vs GL1 spacing;ADC Time Bin [0...1199];Current - previous GL1 GTM BCO [BCO]",
                                          kRun3TruncatedWaveformRecoveryWindowPlotingRange, -.5, static_cast<double>(kRun3TruncatedWaveformRecoveryWindowPlotingRange) - .5, 1001, -.5, 1000.5);
  hm->registerHisto(h_Run3Waveform_GL1Spacing);

  h_Run3WaveformRecovered_GL1Spacing = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3WaveformRecovered_GL1Spacing",  //
                                                   TString(m_HistoPrefix.c_str()) +
                                                       " Run3 matched waveform ADC sum after truncated waveform recovery vs GL1 spacing;ADC Time Bin [0...1199];Current - previous GL1 GTM BCO [BCO]",
                                                   kRun3TruncatedWaveformRecoveryWindowPlotingRange, -.5, static_cast<double>(kRun3TruncatedWaveformRecoveryWindowPlotingRange) - .5, 1001, -.5, 1000.5);
  hm->registerHisto(h_Run3WaveformRecovered_GL1Spacing);

  h_Run3FEE_TimeFrameCount_GL1Spacing = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3FEE_TimeFrameCount_GL1Spacing",  //
                                                 TString(m_HistoPrefix.c_str()) +
                                                     " Run3 exact timeframe count by FEE vs GL1 spacing;Current - previous GL1 GTM BCO [BCO];FEE",
                                                 1001, -.5, 1000.5, MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3FEE_TimeFrameCount_GL1Spacing);

  h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3FEE_TimeFrameRecoveredCount_GL1Spacing",  //
                                                          TString(m_HistoPrefix.c_str()) +
                                                              " Run3 timeframe count by FEE after truncated waveform recovery vs GL1 spacing;Current - previous GL1 GTM BCO [BCO];FEE",
                                                          1001, -.5, 1000.5, MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing);

  h_Run3FEE_TriggerCount_GL1Spacing = new TH2I(TString(m_HistoPrefix.c_str()) + "_Run3FEE_TriggerCount_GL1Spacing",  //
                                               TString(m_HistoPrefix.c_str()) +
                                                   " Run3 GL1 trigger count by FEE vs GL1 spacing;Current - previous GL1 GTM BCO [BCO];FEE",
                                               1001, -.5, 1000.5, MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5);
  hm->registerHisto(h_Run3FEE_TriggerCount_GL1Spacing);

  h_Run3PreviousTimeFrameWaveformADC = new TH1I(TString(m_HistoPrefix.c_str()) + "_Run3PreviousTimeFrameWaveformADCCache",  //
                                                TString(m_HistoPrefix.c_str()) +
                                                    " Run3 previous matched waveform ADC cache;ADC Time Bin [0...1199];Sum ADC",
                                                kRun3TruncatedWaveformRecoveryWindowPlotingRange, -.5, static_cast<double>(kRun3TruncatedWaveformRecoveryWindowPlotingRange) - .5);
  h_Run3PreviousTimeFrameWaveformADC->SetDirectory(nullptr);

  h_Run3PreviousTimeFrameRecoveredWaveformADC = new TH1I(TString(m_HistoPrefix.c_str()) + "_Run3PreviousTimeFrameRecoveredWaveformADCCache",  //
                                                         TString(m_HistoPrefix.c_str()) +
                                                             " Run3 previous matched waveform ADC cache after truncated waveform recovery;ADC Time Bin [0...1199];Sum ADC",
                                                         kRun3TruncatedWaveformRecoveryWindowPlotingRange, -.5, static_cast<double>(kRun3TruncatedWaveformRecoveryWindowPlotingRange) - .5);
  h_Run3PreviousTimeFrameRecoveredWaveformADC->SetDirectory(nullptr);

  h_ProcessPacket_Time = new TH2I(TString(m_HistoPrefix.c_str()) + "_ProcessPacket_Time",  //
                                  TString(m_HistoPrefix.c_str()) +
                                      " Time cost to run ProcessPacket();Call counts;Time elapsed per call [ms];Count",
                                  100, 0, 30e6, 100, 0, 10);
  hm->registerHisto(h_ProcessPacket_Time);
}

TpcTimeFrameBuilderRun3::~TpcTimeFrameBuilderRun3()
{
  for (auto& feeTimeHitMap : m_timeHitMap)
  {
    for (auto& timeHitEntry : feeTimeHitMap)
    {
      for (TpcRawHit* hit : timeHitEntry.second)
      {
        delete hit;
      }
      timeHitEntry.second.clear();
    }
  }

  for (auto& timeFrameEntry : m_timeFrameMap)
  {
    while (!timeFrameEntry.second.empty())
    {
      TpcRawHit* hit = timeFrameEntry.second.back();
      delete hit;
      timeFrameEntry.second.pop_back();
    }
  }

  write_bx_counter_sync_cdb_tree();

  delete h_Run3PreviousTimeFrameWaveformADC;
  delete h_Run3PreviousTimeFrameRecoveredWaveformADC;

  delete m_packetTimer;

  delete m_digitalCurrentDebugTTree;
}

void TpcTimeFrameBuilderRun3::setVerbosity(const int i)
{
  m_verbosity = i;

  for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    bcoMatchingInformation.set_verbosity(i);
  }
}

void TpcTimeFrameBuilderRun3::write_bx_counter_sync_cdb_tree() const
{
  if (m_bxCounterSyncCDBTTreeName.empty())
  {
    return;
  }

  CDBTTree cdbtree(m_bxCounterSyncCDBTTreeName);

  int entry_count = 0;
  for (size_t fee = 0; fee < m_bcoMatchingInformation_vec.size(); ++fee)
  {
    const BcoMatchingInformation& bco_info = m_bcoMatchingInformation_vec[fee];
    const size_t observation_count = std::min(bco_info.get_bx_counter_sync_observation_count(),
                                              BcoMatchingInformation::kMaxBXCounterSyncObservations);
    const auto& observations = bco_info.get_bx_counter_sync_observations();
    for (size_t observation_index = 0; observation_index < observation_count; ++observation_index)
    {
      const BcoMatchingInformation::BXCounterSyncObservation& observation = observations[observation_index];
      const int channel = static_cast<int>(fee * BcoMatchingInformation::kMaxBXCounterSyncObservations + observation_index);
      cdbtree.SetIntValue(channel, "packet_id", m_packet_id);
      cdbtree.SetIntValue(channel, "fee", static_cast<int>(fee));
      cdbtree.SetIntValue(channel, "observation", static_cast<int>(observation_index));
      cdbtree.SetUInt64Value(channel, "bx_counter_sync_gtm_bco", observation.bx_counter_sync_gtm_bco);
      cdbtree.SetUInt64Value(channel, "bco_reference_gtm_bco", observation.bco_reference_gtm_bco);
      cdbtree.SetUInt64Value(channel, "m_bco_reference_gtm_bco", observation.m_bco_reference.first);
      cdbtree.SetIntValue(channel, "m_bco_reference_fee_bco", static_cast<int>(observation.m_bco_reference.second));
      ++entry_count;
    }
  }

  if (entry_count == 0)
  {
    return;
  }

  cdbtree.SetSingleIntValue("packet_id", m_packet_id);
  cdbtree.SetSingleIntValue("n_bx_counter_sync_observations", entry_count);
  cdbtree.SetSingleIntValue("max_fee_count", MAX_FEECOUNT);
  cdbtree.SetSingleIntValue("max_observations_per_fee", static_cast<int>(BcoMatchingInformation::kMaxBXCounterSyncObservations));
  cdbtree.CommitSingle();
  cdbtree.Commit();
  cdbtree.WriteCDBTTree();

  if (m_verbosity >= 0)
  {
    std::cout << __PRETTY_FUNCTION__ << " - saved " << entry_count
              << " BX_COUNTER_SYNC_T observations to " << m_bxCounterSyncCDBTTreeName << std::endl;
  }
}

void TpcTimeFrameBuilderRun3::fill_waveform_gl1_spacing(TH1 *waveform_adc_cache, TH2 *waveform_gl1_spacing, uint64_t gtm_bco_spacing) const
{
  assert(waveform_adc_cache);
  assert(waveform_gl1_spacing);

  const int waveform_ybin = waveform_gl1_spacing->GetYaxis()->FindFixBin(static_cast<double>(gtm_bco_spacing));
  double waveform_entries = 0;
  for (int xbin = 1; xbin <= waveform_adc_cache->GetNbinsX(); ++xbin)
  {
    const double adc_sum = waveform_adc_cache->GetBinContent(xbin);
    if (adc_sum == 0)
    {
      continue;
    }

    waveform_gl1_spacing->AddBinContent(waveform_gl1_spacing->GetBin(xbin, waveform_ybin), adc_sum);
    ++waveform_entries;
  }
  waveform_gl1_spacing->SetEntries(waveform_gl1_spacing->GetEntries() + waveform_entries);
}

void TpcTimeFrameBuilderRun3::flush_previous_timeframe_qa_cache(uint64_t current_gtm_bco)
{
  assert(h_Run3PreviousTimeFrameWaveformADC);
  assert(h_Run3PreviousTimeFrameRecoveredWaveformADC);
  assert(h_Run3Waveform_GL1Spacing);
  assert(h_Run3WaveformRecovered_GL1Spacing);
  assert(h_Run3FEE_TimeFrameCount_GL1Spacing);
  assert(h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing);
  assert(h_Run3FEE_TriggerCount_GL1Spacing);

  if (!m_previousTimeFrameGtmBco)
  {
    return;
  }

  const uint64_t previous_gtm_bco = *m_previousTimeFrameGtmBco;
  static constexpr uint64_t gtm_clock_range = uint64_t(1) << 40U;
  const uint64_t current_gtm_bco_rollover_corrected = current_gtm_bco >= previous_gtm_bco
                                                          ? current_gtm_bco
                                                          : current_gtm_bco + gtm_clock_range;
  const uint64_t gtm_bco_spacing = current_gtm_bco_rollover_corrected - previous_gtm_bco;

  fill_waveform_gl1_spacing(h_Run3PreviousTimeFrameWaveformADC, h_Run3Waveform_GL1Spacing, gtm_bco_spacing);
  fill_waveform_gl1_spacing(h_Run3PreviousTimeFrameRecoveredWaveformADC, h_Run3WaveformRecovered_GL1Spacing, gtm_bco_spacing);

  const int fee_xbin = h_Run3FEE_TriggerCount_GL1Spacing->GetXaxis()->FindFixBin(static_cast<double>(gtm_bco_spacing));
  double timeframe_entries = 0;
  double recovered_timeframe_entries = 0;
  for (uint16_t fee = 0; fee < MAX_FEECOUNT; ++fee)
  {
    const int fee_ybin = static_cast<int>(fee) + 1;
    h_Run3FEE_TriggerCount_GL1Spacing->AddBinContent(h_Run3FEE_TriggerCount_GL1Spacing->GetBin(fee_xbin, fee_ybin));

    if (m_previousTimeFrameExactFees.test(fee))
    {
      h_Run3FEE_TimeFrameCount_GL1Spacing->AddBinContent(h_Run3FEE_TimeFrameCount_GL1Spacing->GetBin(fee_xbin, fee_ybin));
      ++timeframe_entries;
    }

    if (m_previousTimeFrameRecoveredFees.test(fee))
    {
      h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing->AddBinContent(h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing->GetBin(fee_xbin, fee_ybin));
      ++recovered_timeframe_entries;
    }
  }
  h_Run3FEE_TriggerCount_GL1Spacing->SetEntries(h_Run3FEE_TriggerCount_GL1Spacing->GetEntries() + MAX_FEECOUNT);
  h_Run3FEE_TimeFrameCount_GL1Spacing->SetEntries(h_Run3FEE_TimeFrameCount_GL1Spacing->GetEntries() + timeframe_entries);
  h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing->SetEntries(h_Run3FEE_TimeFrameRecoveredCount_GL1Spacing->GetEntries() + recovered_timeframe_entries);

  h_Run3PreviousTimeFrameWaveformADC->Reset();
  h_Run3PreviousTimeFrameRecoveredWaveformADC->Reset();
  m_previousTimeFrameExactFees.reset();
  m_previousTimeFrameRecoveredFees.reset();
  m_previousTimeFrameGtmBco.reset();
}

void TpcTimeFrameBuilderRun3::cache_waveform_adc(TH1 *waveform_adc_cache, const std::vector<TpcRawHit*>& timeframe) const
{
  assert(waveform_adc_cache);

  waveform_adc_cache->Reset();
  for (const TpcRawHit* hit : timeframe)
  {
    if (!hit)
    {
      continue;
    }

    std::unique_ptr<TpcRawHit::AdcIterator> adc_iter(hit->CreateAdcIterator());
    if (!adc_iter)
    {
      continue;
    }

    for (adc_iter->First(); !adc_iter->IsDone(); adc_iter->Next())
    {
      const uint16_t time_bin = adc_iter->CurrentTimeBin();
      const uint16_t adc = adc_iter->CurrentAdc();
      if (adc == 0 || time_bin >= kRun3TruncatedWaveformRecoveryWindow)
      {
        continue;
      }

      waveform_adc_cache->AddBinContent(static_cast<int>(time_bin) + 1, adc);
    }
  }
}

void TpcTimeFrameBuilderRun3::cache_timeframe_qa(uint64_t gtm_bco, const std::vector<TpcRawHit*>& timeframe, const std::bitset<MAX_FEECOUNT>& exact_matched_fees)
{
  cache_waveform_adc(h_Run3PreviousTimeFrameRecoveredWaveformADC, timeframe);

  m_previousTimeFrameExactFees = exact_matched_fees;
  m_previousTimeFrameRecoveredFees.reset();
  for (const TpcRawHit* hit : timeframe)
  {
    if (!hit || hit->get_fee() >= MAX_FEECOUNT)
    {
      continue;
    }
    m_previousTimeFrameRecoveredFees.set(hit->get_fee());
  }
  m_previousTimeFrameGtmBco = gtm_bco;
}

int64_t TpcTimeFrameBuilderRun3::get_signed_fee_bco_diff(uint32_t first, uint32_t second)
{
  static constexpr int64_t fee_clock_range = static_cast<int64_t>(uint64_t{1} << 20U);
  static constexpr int64_t fee_clock_half_range = static_cast<int64_t>(uint64_t{1} << 19U);

  int64_t diff = static_cast<int64_t>(first & kFEEClockMask) - static_cast<int64_t>(second & kFEEClockMask);
  if (diff > fee_clock_half_range)
  {
    diff -= fee_clock_range;
  }
  else if (diff < -fee_clock_half_range)
  {
    diff += fee_clock_range;
  }
  return diff;
}

uint32_t TpcTimeFrameBuilderRun3::get_fee_bco_diff(uint32_t first, uint32_t second)
{
  const int64_t diff = get_signed_fee_bco_diff(first, second);
  return static_cast<uint32_t>(diff < 0 ? -diff : diff);
}

size_t TpcTimeFrameBuilderRun3::move_time_hits(uint32_t fee_bco, uint16_t fee, std::vector<TpcRawHit*>& timeframe)
{
  if (fee >= m_timeHitMap.size())
  {
    return 0;
  }

  auto& fee_time_hits = m_timeHitMap[fee];
  auto it = fee_time_hits.find(fee_bco & kFEEClockMask);
  if (it == fee_time_hits.end())
  {
    return 0;
  }

  std::vector<TpcRawHit*>& hits = it->second;
  const size_t moved = hits.size();
  if (moved == 0)
  {
    fee_time_hits.erase(it);
    return 0;
  }

  timeframe.reserve(timeframe.size() + moved);
  timeframe.insert(timeframe.end(), hits.begin(), hits.end());
  fee_time_hits.erase(it);
  return moved;
}

size_t TpcTimeFrameBuilderRun3::count_time_hits(uint32_t fee_bco, uint16_t fee) const
{
  if (fee >= m_timeHitMap.size())
  {
    return 0;
  }

  const auto& fee_time_hits = m_timeHitMap[fee];
  auto it = fee_time_hits.find(fee_bco & kFEEClockMask);
  if (it == fee_time_hits.end())
  {
    return 0;
  }

  return it->second.size();
}

size_t TpcTimeFrameBuilderRun3::time_hit_bucket_count() const
{
  size_t count = 0;
  for (const auto& fee_time_hits : m_timeHitMap)
  {
    count += fee_time_hits.size();
  }
  return count;
}

std::optional<uint32_t> TpcTimeFrameBuilderRun3::find_fuzzy_fee_bco(uint32_t predicted_fee_bco, uint16_t fee) const
{
  if (fee >= m_timeHitMap.size())
  {
    return std::nullopt;
  }

  const auto& fee_time_hits = m_timeHitMap[fee];
  if (fee_time_hits.empty())
  {
    return std::nullopt;
  }

  predicted_fee_bco &= kFEEClockMask;
  uint32_t best_fee_bco = 0;
  uint32_t best_diff = std::numeric_limits<uint32_t>::max();
  bool found = false;

  auto consider_fee_bco = [&](uint32_t fee_bco, const std::vector<TpcRawHit *>& hits)
  {
    if (hits.empty())
    {
      return;
    }

    const uint32_t diff = get_fee_bco_diff(fee_bco, predicted_fee_bco);
    if (diff <= kRun3FeeMatchWindow && (!found || diff < best_diff || (diff == best_diff && fee_bco < best_fee_bco)))
    {
      found = true;
      best_diff = diff;
      best_fee_bco = fee_bco;
    }
  };

  auto scan_range = [&](uint32_t first_fee_bco, uint32_t last_fee_bco)
  {
    for (auto it = fee_time_hits.lower_bound(first_fee_bco); it != fee_time_hits.end() && it->first <= last_fee_bco; ++it)
    {
      consider_fee_bco(it->first, it->second);
    }
  };

  const uint32_t lower_fee_bco = (predicted_fee_bco - kRun3FeeMatchWindow) & kFEEClockMask;
  const uint32_t upper_fee_bco = (predicted_fee_bco + kRun3FeeMatchWindow) & kFEEClockMask;
  if (lower_fee_bco <= upper_fee_bco)
  {
    scan_range(lower_fee_bco, upper_fee_bco);
  }
  else
  {
    scan_range(lower_fee_bco, kFEEClockMask);
    scan_range(0, upper_fee_bco);
  }

  if (found)
  {
    return best_fee_bco;
  }
  return std::nullopt;
}

size_t TpcTimeFrameBuilderRun3::append_shifted_waveforms(TpcRawHitRun3_typ *target, const TpcRawHit &source, uint32_t fee_clock_shift) const
{
  if (!target || fee_clock_shift >= kRun3TruncatedWaveformRecoveryWindow)
  {
    return 0;
  }

  const auto& source_run3 = static_cast<const TpcRawHitRun3_typ&>(source);
  const TpcRawHitRun3_typ::AdcWaveformVector_t& source_waveforms = source_run3.get_adc_waveforms();
  if (source_waveforms.empty())
  {
    return 0;
  }

  size_t appended_waveforms = 0;
  std::vector<uint16_t> adc_values;
  uint16_t waveform_start = 0;
  uint32_t expected_time_bin = std::numeric_limits<uint32_t>::max();

  auto flush_waveform = [&]()
  {
    if (adc_values.empty())
    {
      return;
    }

    std::vector<uint16_t> waveform_adc;
    waveform_adc.swap(adc_values);
    target->move_adc_waveform(waveform_start, std::move(waveform_adc));
    expected_time_bin = std::numeric_limits<uint32_t>::max();
    ++appended_waveforms;
  };

  for (const TpcRawHitRun3_typ::AdcWaveform_t& source_waveform : source_waveforms)
  {
    const std::vector<uint16_t>& source_adc_values = source_waveform.second;
    if (source_adc_values.empty())
    {
      continue;
    }

    const uint32_t shifted_waveform_start = static_cast<uint32_t>(source_waveform.first) + fee_clock_shift;
    if (shifted_waveform_start >= kRun3TruncatedWaveformRecoveryWindow)
    {
      flush_waveform();
      continue;
    }

    for (size_t adc_index = 0; adc_index < source_adc_values.size(); ++adc_index)
    {
      const uint32_t shifted_time_bin = shifted_waveform_start + static_cast<uint32_t>(adc_index);
      if (shifted_time_bin >= kRun3TruncatedWaveformRecoveryWindow)
      {
        flush_waveform();
        break;
      }

      if (adc_values.empty() || shifted_time_bin != expected_time_bin)
      {
        flush_waveform();
        waveform_start = static_cast<uint16_t>(shifted_time_bin);
      }

      adc_values.push_back(source_adc_values[adc_index]);
      expected_time_bin = shifted_time_bin + 1U;
    }
  }
  flush_waveform();

  return appended_waveforms;
}

size_t TpcTimeFrameBuilderRun3::recover_truncated_waveforms(uint32_t predicted_fee_bco, uint16_t fee, std::vector<TpcRawHit*>& timeframe)
{
  if (fee >= m_timeHitMap.size())
  {
    return 0;
  }

  predicted_fee_bco &= kFEEClockMask;
  std::array<TpcRawHitRun3_typ*, MAX_CHANNELS> current_hits{};
  std::array<int64_t, MAX_CHANNELS> current_hit_diff{};
  current_hit_diff.fill(std::numeric_limits<int64_t>::max());

  for (TpcRawHit* hit : timeframe)
  {
    if (!hit || hit->get_fee() != fee || hit->get_channel() >= MAX_CHANNELS)
    {
      continue;
    }

    TpcRawHitRun3_typ* hit_v3 = dynamic_cast<TpcRawHitRun3_typ*>(hit);
    if (!hit_v3)
    {
      continue;
    }

    const uint16_t channel = hit_v3->get_channel();
    const int64_t signed_diff = get_signed_fee_bco_diff(static_cast<uint32_t>(hit_v3->get_bco()), predicted_fee_bco);
    const int64_t abs_diff = signed_diff < 0 ? -signed_diff : signed_diff;
    if (abs_diff < current_hit_diff[channel])
    {
      current_hit_diff[channel] = abs_diff;
      current_hits[channel] = hit_v3;
    }
  }

  size_t recovered_hits = 0;
  auto& fee_time_hits = m_timeHitMap[fee];
  if (fee_time_hits.empty())
  {
    return 0;
  }

  auto recover_from_bucket = [&](const std::pair<const uint32_t, std::vector<TpcRawHit*>>& bucket)
  {
    for (const TpcRawHit* source_hit : bucket.second)
    {
      if (!source_hit || source_hit->get_channel() >= MAX_CHANNELS)
      {
        continue;
      }

      const uint16_t channel = source_hit->get_channel();
      TpcRawHitRun3_typ* target_hit = current_hits[channel];
      const uint32_t target_fee_bco = target_hit ? static_cast<uint32_t>(target_hit->get_bco()) & kFEEClockMask : predicted_fee_bco;
      const int64_t fee_bco_diff = get_signed_fee_bco_diff(static_cast<uint32_t>(source_hit->get_bco()), target_fee_bco);
      if (fee_bco_diff <= 0 || fee_bco_diff >= static_cast<int64_t>(kRun3TruncatedWaveformRecoveryFEEWindow))
      {
        continue;
      }

      const int64_t fee_clock_shift = fee_bco_diff / static_cast<int64_t>(kRun3FEEClockPerADCClock);
      // FEE BCO is 2x ADC clock, so the waveform shift is half the FEE BCO difference.
      if (fee_clock_shift <= 0)
      {
        continue;
      }

      bool created_target = false;
      if (!target_hit)
      {
        target_hit = new TpcRawHitRun3_typ();
        target_hit->set_bco(predicted_fee_bco);
        target_hit->set_packetid(m_packet_id);
        target_hit->set_fee(fee);
        target_hit->set_channel(channel);
        target_hit->set_type(source_hit->get_type());
        target_hit->set_checksumerror(source_hit->get_checksumerror());
        target_hit->set_parityerror(source_hit->get_parityerror());
        timeframe.push_back(target_hit);
        current_hits[channel] = target_hit;
        current_hit_diff[channel] = 0;
        created_target = true;
      }

      const size_t appended_waveforms = append_shifted_waveforms(target_hit, *source_hit, static_cast<uint32_t>(fee_clock_shift));
      if (appended_waveforms == 0)
      {
        if (created_target)
        {
          current_hits[channel] = nullptr;
          current_hit_diff[channel] = std::numeric_limits<int64_t>::max();
          assert(!timeframe.empty() && timeframe.back() == target_hit);
          timeframe.pop_back();
          delete target_hit;
        }
        continue;
      }

      target_hit->set_checksumerror(target_hit->get_checksumerror() || source_hit->get_checksumerror());
      target_hit->set_parityerror(target_hit->get_parityerror() || source_hit->get_parityerror());
      ++recovered_hits;
    }
  };

  const uint32_t lower_fee_bco = (predicted_fee_bco + 1U) & kFEEClockMask;
  const uint32_t upper_fee_bco = (predicted_fee_bco + kRun3TruncatedWaveformRecoveryFEEWindow - 1U) & kFEEClockMask;
  auto scan_range = [&](uint32_t first_fee_bco, uint32_t last_fee_bco)
  {
    for (auto it = fee_time_hits.lower_bound(first_fee_bco); it != fee_time_hits.end() && it->first <= last_fee_bco; ++it)
    {
      recover_from_bucket(*it);
    }
  };

  if (lower_fee_bco <= upper_fee_bco)
  {
    scan_range(lower_fee_bco, upper_fee_bco);
  }
  else
  {
    scan_range(lower_fee_bco, kFEEClockMask);
    scan_range(0, upper_fee_bco);
  }

  return recovered_hits;
}


void TpcTimeFrameBuilderRun3::cleanup_time_hit_map(uint64_t bclk_rollover_corrected, uint32_t fee_clock_window)
{
  assert(m_hFEEDataStream);

  const size_t nfees = std::min(m_timeHitMap.size(), m_bcoMatchingInformation_vec.size());
  for (size_t fee_index = 0; fee_index < nfees; ++fee_index)
  {
    const uint16_t fee = static_cast<uint16_t>(fee_index);
    auto& fee_time_hits = m_timeHitMap[fee_index];

    const std::optional<uint32_t> predicted_fee_bco = m_bcoMatchingInformation_vec[fee_index].get_predicted_fee_bco(bclk_rollover_corrected);
    if (!predicted_fee_bco)
    {
      if (m_verbosity >= 2)
      {
        std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
            << ": WARNING: No predicted FEE BCO for fee index " << fee_index
            << " with bclk_rollover_corrected: 0x" << std::hex << bclk_rollover_corrected << std::dec
            << ". Clearing time hit map for this fee." << std::endl;
      }

      for (auto map_it = fee_time_hits.begin(); map_it != fee_time_hits.end();)
      {
        for (TpcRawHit* hit : map_it->second)
        {
          m_hFEEDataStream->Fill(fee, "HitUnusedBeforeCleanup", 1);
          delete hit;
        }
        map_it = fee_time_hits.erase(map_it);        
      }

      continue;
    }

    for (auto map_it = fee_time_hits.begin(); map_it != fee_time_hits.end();)
    {
      const int64_t diff = get_signed_fee_bco_diff(map_it->first, *predicted_fee_bco);
      if ((fee_clock_window == 0 && diff <= 0) || (fee_clock_window > 0 && diff < -static_cast<int64_t>(fee_clock_window)))
      {
        for (TpcRawHit* hit : map_it->second)
        {
          m_hFEEDataStream->Fill(fee, "HitUnusedBeforeCleanup", 1);
          delete hit;
        }
        map_it = fee_time_hits.erase(map_it);
      }
      else
      {
        ++map_it;
      }
    }
  }
}

bool TpcTimeFrameBuilderRun3::isMoreDataRequired(const uint64_t& gtm_bco) const
{
  for (const BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    // if (not bcoMatchingInformation.is_verified())
    // {
    //   continue;
    // }

    if (bcoMatchingInformation.isMoreDataRequired(gtm_bco))
    {
      return true;
    }
  }

  if (m_verbosity > 1)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
              << ":PASS: All FEEs satisfied for gtm_bco: 0x" << std::hex << gtm_bco << std::dec << ". Return false."
              << std::endl;
  }
  return false;
}

std::vector<TpcRawHit*>& TpcTimeFrameBuilderRun3::getTimeFrame(const uint64_t& gtm_bco)
{
  assert(m_hNorm);
  const uint64_t bclk_rollover_corrected = m_bcoMatchingInformation_vec[0].get_gtm_rollover_correction(gtm_bco);

  for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    bcoMatchingInformation.cleanup(bclk_rollover_corrected);
  }

  cleanup_time_hit_map(bclk_rollover_corrected, kRun3FeeMatchWindow);

  if (auto cached = m_timeFrameMap.find(bclk_rollover_corrected); cached != m_timeFrameMap.end())
  {
    return cached->second;
  }

  flush_previous_timeframe_qa_cache(bclk_rollover_corrected);

  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "	- packet " << m_packet_id
              << ": getTimeFrame for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
              << ": bclk_rollover_corrected: 0x" << std::hex << bclk_rollover_corrected << std::dec
              << std::endl;
  }

  // Track initial buffer usage
  if (m_verbosity >= 2)
  {
    size_t total_time_hits = 0;
    size_t time_hit_map_buckets = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits += bucket.second.size();
      }
    }
    size_t total_gtm_bco_trig = 0;
    size_t total_bco_heartbeat = 0;
    size_t total_gtm_bco_trigger = 0;
    size_t total_bco_matching = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat += bco_info.get_bco_heartbeat_list_size();
      total_gtm_bco_trigger += bco_info.get_gtm_bco_trigger_map_size();
      total_bco_matching += bco_info.get_bco_matching_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [INITIAL] STL buffer usage - m_timeFrameMap: " << m_timeFrameMap.size()
              << " frames, m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_timeHitMap: " << time_hit_map_buckets << " FEE-BCO buckets, "
              << total_time_hits << " total hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig
              << ", bco_heartbeat: " << total_bco_heartbeat
              << ", gtm_trigger_map: " << total_gtm_bco_trigger
              << ", bco_matching: " << total_bco_matching
              << std::endl;
  }

  auto inserted_frame = m_timeFrameMap.emplace(bclk_rollover_corrected, std::vector<TpcRawHit*>{});
  auto frame_it = inserted_frame.first;
  std::vector<TpcRawHit*>& timeframe = frame_it->second;

  size_t exact_hit_count = 0;
  size_t fallback_hit_count = 0;
  std::bitset<MAX_FEECOUNT> exact_matched_fees;
  std::array<uint32_t, MAX_FEECOUNT> predicted_fee_bcos{};
  std::bitset<MAX_FEECOUNT> predicted_fee_bco_available;

  for (size_t fee_index = 0; fee_index < std::min(m_bcoMatchingInformation_vec.size(), static_cast<size_t>(MAX_FEECOUNT)); ++fee_index)
  {
    const uint16_t fee = static_cast<uint16_t>(fee_index);
    const std::optional<uint32_t> predicted_fee_bco = m_bcoMatchingInformation_vec[fee_index].get_predicted_fee_bco(bclk_rollover_corrected);
    if (!predicted_fee_bco)
    {
      continue;
    }
    predicted_fee_bcos[fee] = *predicted_fee_bco;
    predicted_fee_bco_available.set(fee);

    size_t exact_hits = 0;
    for (int32_t fee_clock_offset = -kRun3ExactMatchWindow; fee_clock_offset <= kRun3ExactMatchWindow; ++fee_clock_offset)
    {
      const uint32_t exact_fee_bco = static_cast<uint32_t>(static_cast<uint64_t>(static_cast<int64_t>(*predicted_fee_bco) + fee_clock_offset) & kFEEClockMask);
      const size_t exact_hits_for_bco = move_time_hits(exact_fee_bco, fee, timeframe);
      if (exact_hits_for_bco == 0)
      {
        continue;
      }

      exact_hits += exact_hits_for_bco;
      assert(h_Run3_FEE_GTMMatching_ClockDiff);
      h_Run3_FEE_GTMMatching_ClockDiff->Fill(static_cast<double>(get_signed_fee_bco_diff(exact_fee_bco, *predicted_fee_bco)),
                                             static_cast<double>(fee),
                                             static_cast<double>(exact_hits_for_bco));
    }

    exact_hit_count += exact_hits;
    if (exact_hits > 0)
    {
      assert(h_Run3TimeFrameExactHit_FEE);
      h_Run3TimeFrameExactHit_FEE->Fill(fee, exact_hits);
      exact_matched_fees.set(fee);
      continue;
    }

    const std::optional<uint32_t> fuzzy_fee_bco = find_fuzzy_fee_bco(*predicted_fee_bco, fee);
    if (!fuzzy_fee_bco)
    {
      continue;
    }


    const size_t fuzzy_hits = count_time_hits(*fuzzy_fee_bco, fee);
    if (fuzzy_hits == 0)
    {
      continue;
    }

    fallback_hit_count += fuzzy_hits;
    assert(h_Run3TimeFrameFuzzyHit_FEE);
    h_Run3TimeFrameFuzzyHit_FEE->Fill(fee, fuzzy_hits);
    assert(h_Run3_FEE_GTMMatching_ClockDiff);
    h_Run3_FEE_GTMMatching_ClockDiff->Fill(static_cast<double>(get_signed_fee_bco_diff(*fuzzy_fee_bco, *predicted_fee_bco)),
                                           static_cast<double>(fee),
                                           static_cast<double>(fuzzy_hits));

    if (m_verbosity >= 2)
    {
      std::cout << __PRETTY_FUNCTION__ << "	- packet " << m_packet_id
                << ": Run3 fuzzy FEE-clock fallback for fee " << fee
                << " predicted 0x" << std::hex << *predicted_fee_bco
                << " matched 0x" << *fuzzy_fee_bco << std::dec
                << " diff " << get_signed_fee_bco_diff(*fuzzy_fee_bco, *predicted_fee_bco)
                << " hits " << fuzzy_hits << std::endl;
    }
  }

  if (fallback_hit_count > 0)
  {
    m_hNorm->Fill("Run3_TimeFrame_FuzzyFallback", 1);
    m_hNorm->Fill("Run3_TimeFrame_FuzzyFallback_Hit_Sum", fallback_hit_count);
  }

  cache_waveform_adc(h_Run3PreviousTimeFrameWaveformADC, timeframe);

  // Track buffer usage after exact and fuzzy hit processing
  if (m_verbosity >= 2)
  {
    size_t total_time_hits_post_exact_fuzzy = 0;
    size_t time_hit_map_buckets_post_exact_fuzzy = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets_post_exact_fuzzy += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits_post_exact_fuzzy += bucket.second.size();
      }
    }
    size_t total_gtm_bco_trig_post = 0;
    size_t total_bco_heartbeat_post = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig_post += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat_post += bco_info.get_bco_heartbeat_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [AFTER EXACT/FUZZY] STL buffer usage - exact_hits: " << exact_hit_count
              << ", fuzzy_hits: " << fallback_hit_count
              << ", timeframe size: " << timeframe.size()
              << ", m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_timeHitMap: " << time_hit_map_buckets_post_exact_fuzzy << " buckets, "
              << total_time_hits_post_exact_fuzzy << " hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_post
              << ", bco_heartbeat: " << total_bco_heartbeat_post << "]"
              << std::endl;
  }

  size_t recovered_hit_count = 0;
  for (uint16_t fee = 0; fee < MAX_FEECOUNT; ++fee)
  {
    if (!predicted_fee_bco_available.test(fee))
    {
      continue;
    }

    const size_t recovered_hits = recover_truncated_waveforms(predicted_fee_bcos[fee], fee, timeframe);
    if (recovered_hits == 0)
    {
      continue;
    }

    recovered_hit_count += recovered_hits;
    if (m_hNormTruncatedWaveformRecoveryFeeFirstBin > 0)
    {
      m_hNorm->Fill(static_cast<double>(m_hNormTruncatedWaveformRecoveryFeeFirstBin + fee), static_cast<double>(recovered_hits));
    }
  }

  if (m_verbosity >= 2 && recovered_hit_count > 0)
  {
    std::cout << __PRETTY_FUNCTION__ << "	- packet " << m_packet_id
              << ": Run3 truncated waveform recovery appended " << recovered_hit_count
              << " later hit segments for gtm_bco: 0x" << std::hex << gtm_bco << std::dec << std::endl;
  }

  // Track buffer usage after recovery
  if (m_verbosity >= 2)
  {
    size_t total_time_hits_post_recovery = 0;
    size_t time_hit_map_buckets_post_recovery = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets_post_recovery += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits_post_recovery += bucket.second.size();
      }
    }
    size_t total_gtm_bco_trig_recovery = 0;
    size_t total_bco_heartbeat_recovery = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig_recovery += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat_recovery += bco_info.get_bco_heartbeat_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [AFTER RECOVERY] STL buffer usage - timeframe size: " << timeframe.size()
              << ", recovered_hits: " << recovered_hit_count
              << ", m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_timeHitMap: " << time_hit_map_buckets_post_recovery << " buckets, "
              << total_time_hits_post_recovery << " hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_recovery
              << ", bco_heartbeat: " << total_bco_heartbeat_recovery << "]"
              << std::endl;
  }

  if (timeframe.empty())
  {
    if (m_verbosity >= 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "	- packet " << m_packet_id
                << ":ERROR: Run3 FEE-clock match failed for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
                << " bclk_rollover_corrected 0x" << std::hex << bclk_rollover_corrected << std::dec
                << ". m_timeHitMap size: " << time_hit_bucket_count() << std::endl;
    }

    if (m_verbosity >= 2)
    {
      size_t total_time_hits_empty = 0;
      size_t time_hit_map_buckets_empty = 0;
      for (const auto& fee_time_hits : m_timeHitMap)
      {
        time_hit_map_buckets_empty += fee_time_hits.size();
        for (const auto& bucket : fee_time_hits)
        {
          total_time_hits_empty += bucket.second.size();
        }
      }
      size_t total_gtm_bco_trig_empty = 0;
      size_t total_bco_heartbeat_empty = 0;
      for (const auto& bco_info : m_bcoMatchingInformation_vec)
      {
        total_gtm_bco_trig_empty += bco_info.get_gtm_bco_trig_list_size();
        total_bco_heartbeat_empty += bco_info.get_bco_heartbeat_list_size();
      }
      std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
                << ": [EMPTY-FRAME ERROR] STL buffer usage - m_timeFrameMap: " << m_timeFrameMap.size()
                << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
                << ", m_timeHitMap: " << time_hit_map_buckets_empty << " buckets, "
                << total_time_hits_empty << " hits"
                << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_empty
                << ", bco_heartbeat: " << total_bco_heartbeat_empty << "]"
                << std::endl;
    }

    m_hNorm->Fill("Run3_TimeFrame_MatchFailed", 1);
    m_hNorm->Fill("GTM_TimeFrame_Unmatched", 1);
    cache_timeframe_qa(bclk_rollover_corrected, timeframe, exact_matched_fees);
    m_timeFrameMap.erase(frame_it);
    static std::vector<TpcRawHit*> empty;
    return empty;
  }

  if (exact_hit_count > 0)
  {
    m_hNorm->Fill("Run3_TimeFrame_Exact_Matched", 1);
  }
  m_hNorm->Fill("GTM_TimeFrame_Matched", 1);
  assert(h_TimeFrame_Matched_Size);
  h_TimeFrame_Matched_Size->Fill(timeframe.size());
  m_hNorm->Fill("GTM_TimeFrame_Matched_Hit_Sum", timeframe.size());
  cache_timeframe_qa(bclk_rollover_corrected, timeframe, exact_matched_fees);
  m_UsedTimeFrameSet.push(bclk_rollover_corrected);

  // Track final buffer usage
  if (m_verbosity >= 2)
  {
    size_t total_time_hits_final = 0;
    size_t time_hit_map_buckets_final = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets_final += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits_final += bucket.second.size();
      }
    }
    size_t total_gtm_bco_trig_final = 0;
    size_t total_bco_heartbeat_final = 0;
    size_t total_gtm_bco_trigger_final = 0;
    size_t total_bco_matching_final = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig_final += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat_final += bco_info.get_bco_heartbeat_list_size();
      total_gtm_bco_trigger_final += bco_info.get_gtm_bco_trigger_map_size();
      total_bco_matching_final += bco_info.get_bco_matching_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [FINAL] STL buffer usage - timeframe size: " << timeframe.size()
              << ", exact_hits: " << exact_hit_count
              << ", fuzzy_hits: " << fallback_hit_count
              << ", recovered_hits: " << recovered_hit_count
              << ", m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_timeHitMap: " << time_hit_map_buckets_final << " buckets, "
              << total_time_hits_final << " hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_final
              << ", bco_heartbeat: " << total_bco_heartbeat_final
              << ", gtm_trigger_map: " << total_gtm_bco_trigger_final
              << ", bco_matching: " << total_bco_matching_final
              << std::endl;
  }

  return timeframe;
}

void TpcTimeFrameBuilderRun3::CleanupUsedPackets(const uint64_t& bclk)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "	- packet " << m_packet_id << ": cleaning up bcos < 0x" << std::hex
              << bclk << std::dec 
              << " and m_UsedTimeFrameSet size: " << m_UsedTimeFrameSet.size()
              << std::endl;
  }

  while (!m_UsedTimeFrameSet.empty())
  {
    const uint64_t bco_completed = m_UsedTimeFrameSet.front();
    m_UsedTimeFrameSet.pop();

    auto it = m_timeFrameMap.find(bco_completed);
    if (it != m_timeFrameMap.end())
    {
      while (!it->second.empty())
      {
        TpcRawHit* hit = it->second.back();
        delete hit;
        it->second.pop_back();
      }
      m_timeFrameMap.erase(it);
    }
  }

  const uint64_t bclk_rollover_corrected = m_bcoMatchingInformation_vec[0].get_gtm_rollover_correction(bclk);
  cleanup_time_hit_map(bclk_rollover_corrected, 0);
}

int TpcTimeFrameBuilderRun3::ProcessPacket(Packet* packet)
{
  static size_t call_count = 0;
  ++call_count;

  if (m_verbosity > 1)
  {
    std::cout << "TpcTimeFrameBuilderRun3::ProcessPacket: " << m_packet_id
              << "\t- Entry " << std::endl;
  }

  if (!packet)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- Error : Invalid packet, doing nothing" << std::endl;
    assert(packet);
    return 0;
  }

  const int packet_hit_format = packet->getHitFormat();
  if (packet_hit_format != IDTPCFEEV5 && packet_hit_format != IDTPCFEEV6)
  {
    std::cout << __PRETTY_FUNCTION__
              << "\t- Error : TpcTimeFrameBuilderRun3 only supports packet formats " << IDTPCFEEV5
              << " or " << IDTPCFEEV6
              << " but received packet format " << packet_hit_format
              << ". Aborting run." << std::endl;
    packet->identify();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_hitFormat < 0)
  {
    m_hitFormat = packet_hit_format;
  }
  else if (packet_hit_format != m_hitFormat)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- Error : packet format changed for packet " << m_packet_id
              << " from " << m_hitFormat << " to " << packet_hit_format
              << ". Aborting run." << std::endl;
    packet->identify();
    return Fun4AllReturnCodes::ABORTRUN;
  }
  assert((packet_hit_format == m_hitFormat));


  if (m_packet_id != packet->getIdentifier())
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- Error : mismatched packet with packet ID expectation of " << m_packet_id << ", but received";
    packet->identify();
    assert(m_packet_id == packet->getIdentifier());
    return 0;
  }

  assert(m_packetTimer);
  if ((m_verbosity == 1 && (call_count % 1000) == 0) || (m_verbosity > 1))
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : received packet ";
    packet->identify();

    m_packetTimer->print_stat();
  }
  m_packetTimer->restart();

  // Track initial buffer usage at start of ProcessPacket
  if (m_verbosity >= 2)
  {
    size_t total_time_hits = 0;
    size_t time_hit_map_buckets = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits += bucket.second.size();
      }
    }
    size_t total_fee_data = 0;
    for (const auto& fee_data_deque : m_feeData)
    {
      total_fee_data += fee_data_deque.size();
    }
    size_t total_gtm_bco_trig = 0;
    size_t total_bco_heartbeat = 0;
    size_t total_gtm_bco_trigger = 0;
    size_t total_bco_matching = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat += bco_info.get_bco_heartbeat_list_size();
      total_gtm_bco_trigger += bco_info.get_gtm_bco_trigger_map_size();
      total_bco_matching += bco_info.get_bco_matching_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [INITIAL] STL buffer usage - m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_feeData total: " << total_fee_data
              << ", m_timeHitMap: " << time_hit_map_buckets << " buckets, " << total_time_hits << " hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig
              << ", bco_heartbeat: " << total_bco_heartbeat
              << ", gtm_trigger_map: " << total_gtm_bco_trigger
              << ", bco_matching: " << total_bco_matching
              << std::endl;
  }

  // //remove after testing
  // ;
  // std::cout <<"packet->lValue(0, N_TAGGER) = "<<packet->lValue(0, "N_TAGGER")<<std::endl;
  // std::cout <<"packet->lValue(0, NR_WF) = "<<packet->iValue(0, "NR_WF")<<std::endl;

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  assert(m_hNorm);
  m_hNorm->Fill("Packet", 1);

  int data_length = packet->getDataLength();  // 32bit length
  assert(h_PacketLength);
  h_PacketLength->Fill(data_length);

  int data_padding = packet->getPadding();  // 32bit padding
  assert(h_PacketLength_Padding);
  h_PacketLength_Padding->Fill(data_padding);
  if (data_padding != 0)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : Warning : suspecious padding "
              << data_padding << "\t- in packet " << m_packet_id << ":" << std::endl;
    packet->identify();
    // packet->dump();
  }

  size_t dma_words_buffer = static_cast<size_t>(data_length) * 2 / DAM_DMA_WORD_LENGTH + 1;
  std::vector<dma_word> buffer(dma_words_buffer);

  int l2 = 0;
  packet->fillIntArray(reinterpret_cast<int*>(buffer.data()), data_length + DAM_DMA_WORD_LENGTH / 2, &l2, "DATA");

  if (data_padding != 0)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- :  data_length = " << data_length
              << "\t- data_padding = " << data_padding << "\t l2 = " << l2 << "\t- in packet " << m_packet_id << ":" << std::endl;
  }

  assert(l2 <= data_length);

  if (l2 < data_padding)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : Error : l2 from fillIntArray() is smaller than padding suggesting an invalid data: " << l2
              << "\t- in packet " << m_packet_id << ". Data length: " << data_length
              << ", data padding: " << data_padding << ". Ignore this packet: " << std::endl;
    packet->identify();
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  l2 -= data_padding;

  assert(l2 >= 0);

  size_t dma_words = static_cast<size_t>(l2) * 2 / DAM_DMA_WORD_LENGTH;
  size_t dma_residual = (static_cast<size_t>(l2) * 2) % DAM_DMA_WORD_LENGTH;
  assert(dma_words <= buffer.size());
  assert(h_PacketLength_Residual);
  h_PacketLength_Residual->Fill(dma_residual);
  if (dma_residual > 0)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : Warning : mismatch of RCDAQ data to DMA transfer. Dropping mismatched data: "
              << dma_residual << "\t- in packet " << m_packet_id << ". Dropping residual data : " << std::endl;

    assert(dma_words + 1 < buffer.size());
    const dma_word& last_dma_word_data = buffer[dma_words + 1];
    const uint16_t* last_dma_word = reinterpret_cast<const uint16_t*>(&last_dma_word_data);

    for (size_t i = 0; i < dma_residual; ++i)
    {
      std::cout << "\t- 0x" << std::hex << last_dma_word[i] << std::dec;
    }
    std::cout << std::endl;
  }

  if (m_verbosity > 1)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : packet" << m_packet_id << std::endl
              << "\t-   data_length = " << data_length << std::endl
              << "\t-   data_padding = " << data_padding << std::endl
              << "\t-   dma_words_buffer = " << dma_words_buffer << std::endl
              << "\t-   l2 = " << l2 << std::endl
              << "\t-   dma_words = " << dma_words << std::endl;
  }

  // demultiplexer
  for (size_t index = 0; index < dma_words; ++index)
  {
    const dma_word& dma_word_data = buffer[index];

    if (m_verbosity > 2)
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- : processing DMA word "
                << index << "/" << dma_words << "\t- with header 0x"
                << std::hex << dma_word_data.dma_header << std::dec << std::endl;
    }

    if ((dma_word_data.dma_header & 0xFF00U) == FEE_MAGIC_KEY)
    {
      unsigned int fee_id = dma_word_data.dma_header & 0xffU;

      // for packet id 4XYZ ebdc is XY, endpoint is Z
      if (m_maskedFEEs[((m_packet_id / 10) % 100)].contains(fee_id))
      {
        continue;
      }

      if (fee_id < MAX_FEECOUNT)
      {
        for (const uint16_t& i : dma_word_data.data)
        {
          m_feeData[fee_id].push_back(i);
        }
        m_hNorm->Fill("DMA_WORD_FEE", 1);

        // immediate fee buffer processing to reduce memory consuption
        process_fee_data(fee_id);
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE ID " << fee_id << "\t- at position " << index << std::endl;
        index += DAM_DMA_WORD_LENGTH - 1;
        m_hNorm->Fill("DMA_WORD_FEE_INVALID", 1);
      }
    }

    else if ((dma_word_data.dma_header & 0xFF00U) == GTM_MAGIC_KEY)
    {
      decode_gtm_data(dma_word_data);
      m_hNorm->Fill("DMA_WORD_GTM", 1);
    }
    else
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Unknown data type at position " << index << ": "
                << std::hex << buffer[index].dma_header << std::dec << std::endl;
      // not FEE data, e.g. GTM data or other stream, to be decoded
      m_hNorm->Fill("DMA_WORD_INVALID", 1);
    }
  }

  // Track buffer usage after DMA word processing
  if (m_verbosity >= 2)
  {
    size_t total_time_hits_post = 0;
    size_t time_hit_map_buckets_post = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets_post += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits_post += bucket.second.size();
      }
    }
    size_t total_fee_data_post = 0;
    for (const auto& fee_data_deque : m_feeData)
    {
      total_fee_data_post += fee_data_deque.size();
    }
    size_t total_gtm_bco_trig_post = 0;
    size_t total_bco_heartbeat_post = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig_post += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat_post += bco_info.get_bco_heartbeat_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [AFTER DMA PROCESSING] STL buffer usage - m_feeData total: " << total_fee_data_post
              << ", m_timeHitMap: " << time_hit_map_buckets_post << " buckets, " << total_time_hits_post << " hits"
              << ", m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_post
              << ", bco_heartbeat: " << total_bco_heartbeat_post << "]"
              << std::endl;
  }

  // sanity check for the cached FEE-clock hit size
  for (size_t fee = 0; fee < m_timeHitMap.size(); ++fee)
  {
    auto& fee_time_hits = m_timeHitMap[fee];
    for (auto timehit = fee_time_hits.begin(); timehit != fee_time_hits.end();)
    {
      if (timehit->second.size() > kMaxRawHitLimit)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : Warning : impossible amount of hits for FEE "
                  << fee << " at FEE BCO " << timehit->first << "\t- : " << timehit->second.size()
                  << ", limit is " << kMaxRawHitLimit
                  << ". Dropping this FEE-clock cache!"
                  << std::endl;
        m_hNorm->Fill("TimeFrameSizeLimitError", 1);

        for (TpcRawHit* hit : timehit->second)
        {
          delete hit;
        }
        timehit = fee_time_hits.erase(timehit);
      }
      else
      {
        ++timehit;
      }
    }
  }

  m_packetTimer->stop();
  assert(h_ProcessPacket_Time);
  h_ProcessPacket_Time->Fill(call_count, m_packetTimer->elapsed());

  // Track final buffer usage at end of ProcessPacket
  if (m_verbosity >= 1)
  {
    size_t total_time_hits_final = 0;
    size_t time_hit_map_buckets_final = 0;
    for (const auto& fee_time_hits : m_timeHitMap)
    {
      time_hit_map_buckets_final += fee_time_hits.size();
      for (const auto& bucket : fee_time_hits)
      {
        total_time_hits_final += bucket.second.size();
      }
    }
    size_t total_fee_data_final = 0;
    for (const auto& fee_data_deque : m_feeData)
    {
      total_fee_data_final += fee_data_deque.size();
    }
    size_t total_gtm_bco_trig_final = 0;
    size_t total_bco_heartbeat_final = 0;
    size_t total_gtm_bco_trigger_final = 0;
    size_t total_bco_matching_final = 0;
    for (const auto& bco_info : m_bcoMatchingInformation_vec)
    {
      total_gtm_bco_trig_final += bco_info.get_gtm_bco_trig_list_size();
      total_bco_heartbeat_final += bco_info.get_bco_heartbeat_list_size();
      total_gtm_bco_trigger_final += bco_info.get_gtm_bco_trigger_map_size();
      total_bco_matching_final += bco_info.get_bco_matching_list_size();
    }
    std::cout << __PRETTY_FUNCTION__ << " - packet " << m_packet_id
              << ": [FINAL] STL buffer usage - m_timeFrameMap: " << m_timeFrameMap.size()
              << ", m_UsedTimeFrameSet: " << m_UsedTimeFrameSet.size()
              << ", m_feeData total: " << total_fee_data_final
              << ", m_timeHitMap: " << time_hit_map_buckets_final << " buckets, " << total_time_hits_final << " hits"
              << ", BcoMatchingInfo[gtm_trig: " << total_gtm_bco_trig_final
              << ", bco_heartbeat: " << total_bco_heartbeat_final
              << ", gtm_trigger_map: " << total_gtm_bco_trigger_final
              << ", bco_matching: " << total_bco_matching_final
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTimeFrameBuilderRun3::process_fee_data(unsigned int fee)
{
  assert(m_hFEEDataStream);

  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : processing FEE " << fee << "\t- with " << m_feeData[fee].size() << "\t- words" << std::endl;
  }

  assert(fee < m_feeData.size());
  std::deque<uint16_t>& data_buffer = m_feeData[fee];

  while (HEADER_LENGTH <= data_buffer.size())
  {
    // packet loop

    bool is_digital_current = false;
    // test if digital current packet
    if (data_buffer[3] == FEE_PACKET_MAGIC_KEY_3_DC)
    {
      if (m_verbosity > 2)
      {
        std::cout << __PRETTY_FUNCTION__
                  << "\t- : processing FEE " << fee
                  << "\t- with digital packet" << std::endl;
      }

      m_hFEEDataStream->Fill(fee, "WordDigitalCurrentKeyWord", 1);
      is_digital_current = true;
    }  //     if (data_buffer[3] == FEE_PACKET_MAGIC_KEY_3)
    else
    {
      if (data_buffer[1] != FEE_PACKET_MAGIC_KEY_1)
      {
        if (m_verbosity > 1)
        {
          std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE magic key at position 1 0x" << std::hex << data_buffer[1] << std::dec << std::endl;
        }
        m_hFEEDataStream->Fill(fee, "WordSkipped", 1);
        data_buffer.pop_front();
        continue;
      }
      assert(data_buffer[1] == FEE_PACKET_MAGIC_KEY_1);

      if (data_buffer[2] != FEE_PACKET_MAGIC_KEY_2)
      {
        if (m_verbosity > 1)
        {
          std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE magic key at position 2 0x" << std::hex << data_buffer[2] << std::dec << std::endl;
        }
        m_hFEEDataStream->Fill(fee, "WordSkipped", 1);
        data_buffer.pop_front();
        continue;
      }
      assert(data_buffer[2] == FEE_PACKET_MAGIC_KEY_2);
    }

    // valid packet
    const uint16_t pkt_length = data_buffer[0];  // this is indeed the number of 10-bit words + 5 in this packet
    if (pkt_length > MAX_PACKET_LENGTH)
    {
      if (m_verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE pkt_length " << pkt_length << std::endl;
      }
      m_hFEEDataStream->Fill(fee, "InvalidLength", 1);
      data_buffer.pop_front();
      continue;
    }

    if (pkt_length + 1U > data_buffer.size())
    {
      if (m_verbosity > 2)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : packet over buffer boundary for now, skip decoding and wait for more data: "
                                            " pkt_length = "
                  << pkt_length
                  << "\t- data_buffer.size() = " << data_buffer.size()
                  << std::endl;
      }
      break;
    }

    if (is_digital_current)
    {
      process_fee_data_digital_current(fee, data_buffer);
    }
    else
    {
      process_fee_data_waveform(fee, data_buffer);
    }
    data_buffer.erase(data_buffer.begin(), data_buffer.begin() + pkt_length + 1);
    m_hFEEDataStream->Fill(fee, "WordValid", pkt_length + 1);

  }  //     while (HEADER_LENGTH < data_buffer.size())

  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcTimeFrameBuilderRun3::process_fee_data_waveform(const unsigned int& fee, std::deque<uint16_t>& data_buffer)
{
  const uint16_t& pkt_length = data_buffer[0];

  fee_payload payload;
  // continue the decoding
  payload.fee_id = fee;
  payload.adc_length = data_buffer[0] - HEADER_LENGTH;  // this is indeed the number of 10-bit words in this packet
  payload.data_parity = data_buffer[4] >> 9U;
  payload.sampa_address = static_cast<uint16_t>(data_buffer[4] >> 5U) & 0xfU;
  payload.sampa_channel = data_buffer[4] & 0x1fU;
  payload.channel = data_buffer[4] & 0x1ffU;
  payload.type = static_cast<uint16_t>(data_buffer[3] >> 7U) & 0x7U;
  payload.user_word = data_buffer[3] & 0x7fU;
  payload.bx_timestamp = static_cast<uint32_t>(static_cast<uint32_t>(data_buffer[6] & 0x3ffU) << 10U) | (data_buffer[5] & 0x3ffU);
  payload.data_crc = data_buffer[pkt_length];

  if (!m_fastBCOSkip)
  {
    auto crc_parity = crc16_parity(fee, pkt_length);
    payload.calc_crc = crc_parity.first;
    payload.calc_parity = crc_parity.second;

    if (payload.data_crc != payload.calc_crc)
    {
      if (m_verbosity > 2)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : CRC error in FEE "
                  << fee << "\t- at position " << pkt_length - 1
                  << ": data_crc = " << payload.data_crc
                  << "\t- calc_crc = " << payload.calc_crc << std::endl;
      }
      m_hFEEDataStream->Fill(fee, "HitCRCError", 1);
      // continue;
    }

    if (payload.data_parity != payload.calc_parity)
    {
      if (m_verbosity > 2)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : parity error in FEE "
                  << fee << "\t- at position " << pkt_length - 1
                  << ": data_parity = " << payload.data_parity
                  << "\t- calc_parity = " << payload.calc_parity << std::endl;
      }
      m_hFEEDataStream->Fill(fee, "ParityError", 1);
      // continue;
    }
  }  //     if (not m_fastBCOSkip)

  assert(fee < m_bcoMatchingInformation_vec.size());
  BcoMatchingInformation& m_bcoMatchingInformation = m_bcoMatchingInformation_vec[fee];
  // gtm_bco matching
  if (payload.type == TpcTimeFrameBuilderRun3::BcoMatchingInformation::HEARTBEAT_T)
  {
    if (m_verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__
                << "\t- : received heartbeat packet from FEE " << fee << std::endl;
    }

    // if bco matching information is still not verified, drop the packet
    if (!m_bcoMatchingInformation.is_verified())
    {
      m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncUnavailable", 1);

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilderRun3::process_fee_data - bco_matching not verified for heart beat, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
    }
    else  //       if (not m_bcoMatchingInformation.is_verified())
    {
      const std::optional<uint64_t> result = m_bcoMatchingInformation.find_reference_heartbeat(payload);
      m_hFEEDataStream->Fill(fee, "PacketHeartBeat", 1);

      if (result)
      {
        // assign gtm bco
        payload.gtm_bco = result.value();
        payload.has_clock_sync = true;
        m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncOK", 1);

        assert(m_hFEESAMPAHeartBeatSync);
        m_hFEESAMPAHeartBeatSync->Fill(fee * MAX_SAMPA + payload.sampa_address, 1);
      }
      else
      {
        m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncError", 1);

        // skip the waverform
      }
      if (m_verbosity > 2)
      {
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
    }
  }
  else if (!m_fastBCOSkip)  //     if (payload.type == m_bcoMatchingInformation.HEARTBEAT_T)
  {
    m_hFEEChannelPacketCount->Fill(fee * MAX_CHANNELS + payload.channel, 1);

    // if bco matching information is still not verified, drop the packet
    if (!m_bcoMatchingInformation.is_verified())
    {
      m_hFEEDataStream->Fill(fee, "PacketClockSyncUnavailable", 1);

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilderRun3::process_fee_data - bco_matching not verified, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
    }
    else
    {
      payload.has_clock_sync = true;
      m_hFEEDataStream->Fill(fee, "PacketClockSyncOK", 1);
    }
  }

  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : received data packet "
              << "\t- from FEE " << fee << std::endl
              << "\t- pkt_length = " << pkt_length << std::endl
              << "\t- type = " << payload.type << std::endl
              << "\t- adc_length = " << payload.adc_length << std::endl
              << "\t- sampa_address = " << payload.sampa_address << std::endl
              << "\t- sampa_channel = " << payload.sampa_channel << std::endl
              << "\t- channel = " << payload.channel << std::endl
              << "\t- bx_timestamp = 0x" << std::hex << payload.bx_timestamp << std::dec << std::endl
              << "\t- bco = 0x" << std::hex << payload.gtm_bco << std::dec << std::endl
              << "\t- data_crc = 0x" << std::hex << payload.data_crc << std::dec << std::endl
              << "\t- calc_crc = 0x" << std::hex << payload.calc_crc << std::dec << std::endl
              << "\t- data_parity = 0x" << std::hex << payload.data_parity << std::dec << std::endl
              << "\t- calc_parity = 0x" << std::hex << payload.calc_parity << std::dec << std::endl;
  }

  if ((!m_fastBCOSkip) && payload.has_clock_sync)
  {
    m_hFEEDataStream->Fill(fee, "RawHit", 1);

    // Format is (N sample) (start time), (1st sample)... (Nth sample)
    size_t pos = HEADER_LENGTH;
    std::deque<uint16_t>::const_iterator data_buffer_iterator = data_buffer.cbegin();
    std::advance(data_buffer_iterator, pos);
    while (pos + 2 < pkt_length)
    {
      const uint16_t& nsamp = *data_buffer_iterator;
      ++pos;
      ++data_buffer_iterator;
      const uint16_t& start_t = *data_buffer_iterator;
      ++pos;
      ++data_buffer_iterator;
      if (m_verbosity > 3)
      {
        std::cout << __PRETTY_FUNCTION__ << ": nsamp: " << nsamp
                  << "+ pos: " << pos
                  << " pkt_length: " << pkt_length << " start_t:" << start_t << std::endl;
      }

      if (pos + nsamp > pkt_length)
      {
        if (m_verbosity > 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": WARNING : nsamp: " << nsamp
                    << "+ pos: " << pos
                    << " > pkt_length: " << pkt_length << ", format error over length: " << std::endl;

          for (int print_pos = 0; print_pos <= pkt_length; ++print_pos)
          {
            std::cout << "\t[" << print_pos << "]=0x" << std::hex << data_buffer[print_pos] << std::dec << "(" << data_buffer[print_pos] << ")";
          }
          std::cout << std::endl;
        }
        m_hFEEDataStream->Fill(fee, "HitFormatErrorOverLength", 1);

        break;
      }

      const unsigned int fee_sampa_address = fee * MAX_SAMPA + payload.sampa_address;
      std::vector<uint16_t> adc(nsamp);
      for (int j = 0; j < nsamp; j++)
      {
        const uint16_t& adc_value = *data_buffer_iterator;

        adc[j] = adc_value;
        m_hFEESAMPAADC->Fill(start_t + j, fee_sampa_address, adc_value);

        ++pos;
        ++data_buffer_iterator;  // data_buffer[pos++];
      }
      payload.waveforms.emplace_back(start_t, std::move(adc));

      //   // an exception to deal with the last sample that is missing in the current hit format
      //   if (pos + 1 == pkt_length) break;
    }

    if (pos != pkt_length)
    {
      if (m_verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << ": WARNING : residual data at the end of decoding:"
                  << " pos: " << pos
                  << " <pkt_length: " << pkt_length << ", format error under length" << std::endl;
      }
      m_hFEEDataStream->Fill(fee, "HitFormatErrorMismatchedLength", 1);
    }

    // valid packet in the buffer, create a new hit
    if (payload.type != TpcTimeFrameBuilderRun3::BcoMatchingInformation::HEARTBEAT_T)
    {
      if (fee >= m_timeHitMap.size())
      {
        std::cout << __PRETTY_FUNCTION__ << ": ERROR : invalid FEE " << fee
                  << " for packet " << m_packet_id << ". Dropping waveform hit." << std::endl;
        return;
      }

      TpcRawHitRun3_typ* hit = new TpcRawHitRun3_typ();

      hit->set_bco(payload.bx_timestamp);
      hit->set_packetid(m_packet_id);
      hit->set_fee(fee);
      hit->set_channel(payload.channel);
      hit->set_type(payload.type);
      // hit->set_checksum(payload.data_crc);
      hit->set_checksumerror(payload.data_crc != payload.calc_crc);
      // hit->set_parity(payload.data_parity);
      hit->set_parityerror(payload.data_parity != payload.calc_parity);

      for (std::pair<uint16_t, std::vector<uint16_t>>& waveform : payload.waveforms)
      {
        hit->move_adc_waveform(waveform.first, std::move(waveform.second));
      }

      m_timeHitMap[fee][payload.bx_timestamp & kFEEClockMask].push_back(hit);
    }
  }  //     if (not m_fastBCOSkip)

  return;
}

void TpcTimeFrameBuilderRun3::process_fee_data_digital_current(const unsigned int& fee, std::deque<uint16_t>& data_buffer)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : processing digital_current data " << std::endl;
  }
  m_hFEEDataStream->Fill(fee, "DigitalCurrent", 1);
  const uint16_t& pkt_length = data_buffer[0];

  if (pkt_length != HEADER_LENGTH + digital_current_payload::MAX_CHANNELS * 2 * 2)
  {
    if (m_verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE pkt_length " << pkt_length
                << ", expected at least " << HEADER_LENGTH + digital_current_payload::MAX_CHANNELS * 2 * 2
                << std::endl;
    }
    m_hFEEDataStream->Fill(fee, "DigitalCurrentFormatErrorMismatchedLength", 1);
    return;
  }

  digital_current_payload payload;

  payload.fee = fee;
  payload.pkt_length = pkt_length;
  payload.sampa_address = (data_buffer[4] >> 5U) & 0xfU;  // NOLINT(hicpp-signed-bitwise)
  // payload.sampa_max_channel = data_buffer[4] & 0x1fU;
  payload.channel = data_buffer[4] & 0x1ffU;
  // payload.type          = data_buffer[3];
  payload.bx_timestamp = ((data_buffer[6] & 0x3ffU) << 10U) | (data_buffer[5] & 0x3ff);  // NOLINT(hicpp-signed-bitwise)

  uint16_t pos = HEADER_LENGTH;
  for (int ich = 0; ich < digital_current_payload::MAX_CHANNELS; ich++)
  {
    payload.current[ich] = ((unsigned int) data_buffer[pos]) << 16U | ((unsigned int) data_buffer[pos + 1U]);
    pos++;
    pos++;
    payload.nsamples[ich] = ((unsigned int) data_buffer[pos]) << 16U | ((unsigned int) data_buffer[pos + 1U]);
    pos++;
    pos++;
  }

  if (pos != pkt_length)
  {
    if (m_verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- : Warning : residual data at the end of decoding:"
                << " pos: " << pos
                << " <pkt_length: " << pkt_length << ", format error under length" << std::endl;
    }
  }

  payload.data_crc = data_buffer[pkt_length];
  auto crc_parity = crc16_parity(fee, pkt_length);
  payload.calc_crc = crc_parity.first;
  // payload.calc_parity = crc_parity.second;

  if (payload.data_crc != payload.calc_crc)
  {
    if (m_verbosity > 2)
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- : CRC error in FEE "
                << fee << "\t- at position " << pkt_length - 1
                << ": data_crc = " << payload.data_crc
                << "\t- calc_crc = " << payload.calc_crc << std::endl;
    }
    m_hFEEDataStream->Fill(fee, "DigitalCurrentCRCError", 1);
    // continue;
  }

  assert(fee < m_bcoMatchingInformation_vec.size());
  BcoMatchingInformation& m_bcoMatchingInformation = m_bcoMatchingInformation_vec[fee];
  std::tie(payload.gtm_bco, payload.bx_timestamp_predicted) = m_bcoMatchingInformation.find_dc_read_bco();

  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : received digital current packet "
              << "\t- from FEE " << fee << std::endl
              << "\t- pkt_length = " << pkt_length << std::endl
              << "\t- sampa_address = " << payload.sampa_address << std::endl
              << "\t- channel = " << payload.channel << std::endl
              << "\t- bx_timestamp = 0x" << std::hex << payload.bx_timestamp << std::dec << std::endl
              << "\t- gtm_bco = 0x" << std::hex << payload.gtm_bco << std::dec << std::endl
              << "\t- bx_timestamp_predicted = 0x" << std::hex << payload.bx_timestamp_predicted << std::dec << std::endl;

    std::cout << "\t- current:";
    for (int ich = 0; ich < digital_current_payload::MAX_CHANNELS; ich++)
    {
      std::cout << "\t[" << ich << "] = " << payload.current[ich];
    }
    std::cout << std::endl;
    std::cout << "\t- nsamples:";
    for (int ich = 0; ich < digital_current_payload::MAX_CHANNELS; ich++)
    {
      std::cout << "\t[" << ich << "] = " << payload.nsamples[ich];
    }
    std::cout << std::endl;
    std::cout << "\t- data_crc = 0x" << std::hex << payload.data_crc << std::dec << std::endl
              << "\t- calc_crc = 0x" << std::hex << payload.calc_crc << std::dec << std::endl;
  }

  if (m_digitalCurrentDebugTTree)
  {
    m_digitalCurrentDebugTTree->fill(payload);
  }

  return;
}

void TpcTimeFrameBuilderRun3::SaveBXCounterSyncCDBTTree(const std::string& name)
{
  m_bxCounterSyncCDBTTreeName = name;

  if (m_verbosity >= 1)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : Saving BX counter sync CDB TTree to " << m_bxCounterSyncCDBTTreeName << std::endl;
  }
}

void TpcTimeFrameBuilderRun3::SaveDigitalCurrentDebugTTree(const std::string& name)
{
  if (m_verbosity >= 1)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : Saving digital current debug TTree to " << name << std::endl;
  }

  m_digitalCurrentDebugTTree = new TpcTimeFrameBuilderRun3::DigitalCurrentDebugTTree(name);
}

TpcTimeFrameBuilderRun3::DigitalCurrentDebugTTree::DigitalCurrentDebugTTree(const std::string& name)
  : m_name(name)
{
  // open TFile
  PHTFileServer::open(m_name, "RECREATE");

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_tDigitalCurrent = new TTree("T_DigitalCurrent", "DigitalCurrent Debug TTree");
  assert(m_tDigitalCurrent);

  m_tDigitalCurrent->Branch("dc", &m_payload,
                            "gtm_bco/l:bx_timestamp_predicted/i:fee/s:pkt_length/s:channel/s:sampa_address/s:bx_timestamp/i:current[8]/i:nsamples[8]/i:data_crc/s:calc_crc/s");
}

TpcTimeFrameBuilderRun3::DigitalCurrentDebugTTree::~DigitalCurrentDebugTTree()
{
  // open TFile
  PHTFileServer::write(m_name);
}

void TpcTimeFrameBuilderRun3::DigitalCurrentDebugTTree::fill(const TpcTimeFrameBuilderRun3::digital_current_payload& payload)
{
  assert(m_tDigitalCurrent);

  m_payload = payload;
  m_tDigitalCurrent->Fill();
}

int TpcTimeFrameBuilderRun3::decode_gtm_data(const TpcTimeFrameBuilderRun3::dma_word& gtm_word)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- : processing GTM data " << std::endl;
  }

  const uint8_t* gtm = reinterpret_cast<const uint8_t*>(&gtm_word);

  gtm_payload payload;

  payload.pkt_type = gtm[0] | static_cast<uint16_t>((unsigned short) gtm[1] << 8U);
  //    if (payload.pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY && payload.pkt_type != GTM_ENDAT_MAGIC_KEY)
  if (payload.pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY && payload.pkt_type != GTM_ENDAT_MAGIC_KEY && payload.pkt_type != GTM_MODEBIT_MAGIC_KEY)
  {
    return -1;
  }

  payload.is_lvl1 = payload.pkt_type == GTM_LVL1_ACCEPT_MAGIC_KEY;
  payload.is_endat = payload.pkt_type == GTM_ENDAT_MAGIC_KEY;
  payload.is_modebit = payload.pkt_type == GTM_MODEBIT_MAGIC_KEY;

  payload.bco = ((unsigned long long) gtm[2] << 0U) | ((unsigned long long) gtm[3] << 8U) | ((unsigned long long) gtm[4] << 16U) | ((unsigned long long) gtm[5] << 24U) | ((unsigned long long) gtm[6] << 32U) | (((unsigned long long) gtm[7]) << 40U);
  payload.lvl1_count = ((unsigned int) gtm[8] << 0U) | ((unsigned int) gtm[9] << 8U) | ((unsigned int) gtm[10] << 16U) | ((unsigned int) gtm[11] << 24U);
  payload.endat_count = ((unsigned int) gtm[12] << 0U) | ((unsigned int) gtm[13] << 8U) | ((unsigned int) gtm[14] << 16U) | ((unsigned int) gtm[15] << 24U);
  payload.last_bco = ((unsigned long long) gtm[16] << 0U) | ((unsigned long long) gtm[17] << 8U) | ((unsigned long long) gtm[18] << 16U) | ((unsigned long long) gtm[19] << 24U) | ((unsigned long long) gtm[20] << 32U) | (((unsigned long long) gtm[21]) << 40U);
  payload.modebits = gtm[22];
  payload.userbits = gtm[23];

  if (m_verbosity >= 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- GTM data : "
              << "\t- pkt_type = " << payload.pkt_type << std::endl
              << "\t- is_lvl1 = " << payload.is_lvl1 << std::endl
              << "\t- is_endat = " << payload.is_endat << std::endl
              << "\t- is_modebit = " << payload.is_modebit << std::endl
              << "\t- bco = 0x" << std::hex << payload.bco << std::dec << std::endl
              << "\t- lvl1_count = " << payload.lvl1_count << std::endl
              << "\t- endat_count = " << payload.endat_count << std::endl
              << "\t- last_bco = 0x" << std::hex << payload.last_bco << std::dec << std::endl
              << "\t- modebits =  0x" << std::hex << (int) payload.modebits << std::dec << std::endl
              << "\t- userbits =  0x" << std::hex << (int) payload.userbits << std::dec << std::endl;
  }

  if (payload.is_modebit)
  {
    if (payload.modebits == BcoMatchingInformation::ELINK_HEARTBEAT_T)
    {
      if (m_verbosity > 2)
      {
        std::cout << "\t- (Heartbeat modebit)" << std::endl;
      }
      assert(m_hNorm);
      m_hNorm->Fill("DMA_WORD_GTM_HEARTBEAT", 1);
    }

    if (payload.modebits == BcoMatchingInformation::DC_STOP_SEND_T)
    {
      if (m_verbosity > 2)
      {
        std::cout << "\t- (DC stop send modebit)" << std::endl;
      }
      assert(m_hNorm);
      m_hNorm->Fill("DMA_WORD_GTM_DC_STOP_SEND", 1);
    }
  }

  if (!(m_fastBCOSkip && (payload.is_lvl1 || payload.is_endat)))
  {
    int fee = -1;
    for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
    {
      ++fee;

      if (m_verbosity > 2)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- : processing GTM data for FEE " << fee << std::endl;
      }

      bcoMatchingInformation.save_gtm_bco_information(payload);

      if (m_verbosity > 2)
      {
        bcoMatchingInformation.print_gtm_bco_information();
      }
    }
  }  //   if (not m_fastBCOSkip)

  return 0;
}

uint16_t TpcTimeFrameBuilderRun3::reverseBits(const uint16_t x) const
{
  uint16_t n = x;
  n = (static_cast<uint16_t>(n >> 1U) & 0x55555555U) | (static_cast<uint16_t>(n << 1U) & 0xaaaaaaaaU);
  n = (static_cast<uint16_t>(n >> 2U) & 0x33333333U) | (static_cast<uint16_t>(n << 2U) & 0xccccccccU);
  n = (static_cast<uint16_t>(n >> 4U) & 0x0f0f0f0fU) | (static_cast<uint16_t>(n << 4U) & 0xf0f0f0f0U);
  n = (static_cast<uint16_t>(n >> 8U) & 0x00ff00ffU) | (static_cast<uint16_t>(n << 8U) & 0xff00ff00U);
  // n = (n >> 16U) & 0x0000ffffU | (n << 16U) & 0xffff0000U;
  return n;
}

std::pair<uint16_t, uint16_t> TpcTimeFrameBuilderRun3::crc16_parity(const uint32_t fee, const uint16_t l) const
{
  const std::deque<uint16_t>& data_buffer = m_feeData[fee];
  assert(l < data_buffer.size());

  std::deque<uint16_t>::const_iterator it = data_buffer.begin();

  uint16_t crc = 0xffffU;
  uint16_t data_parity = 0U;

  for (int i = 0; i < l; ++i, ++it)
  {
    const uint16_t& x = *it;

    crc ^= reverseBits(x);
    for (uint16_t k = 0; k < 16U; k++)
    {
      crc = crc & 1U ? static_cast<uint16_t>(crc >> 1U) ^ 0xa001U : crc >> 1U;
    }

    // parity on data payload only
    if (i >= HEADER_LENGTH)
    {
      // fast parity
      uint16_t word = x & uint16_t((1U << 10U) - 1U);
      word = word ^ static_cast<uint16_t>(word >> 1U);
      word = word ^ static_cast<uint16_t>(word >> 2U);
      word = word ^ static_cast<uint16_t>(word >> 4U);
      word = word ^ static_cast<uint16_t>(word >> 8U);
      data_parity ^= word & 1U;
    }
  }
  crc = reverseBits(crc);
  return std::make_pair(crc, data_parity);
}

namespace
{
  // streamer for lists
  template <class T>
  std::ostream& operator<<(std::ostream& o, const std::list<T>& list)
  {
    if (list.empty())
    {
      o << "{}";
    }
    else
    {
      const bool is_hex = (o.flags() & std::ios_base::hex);
      o << "{ ";
      bool first = true;
      for (const auto& value : list)
      {
        if (!first)
        {
          o << ", ";
        }
        if (is_hex)
        {
          o << "0x";
        }
        o << value;
        first = false;
      }
      o << "\t- }";
    }
    return o;
  }

  template <class T>
  std::ostream& operator<<(std::ostream& o, const std::vector<T>& list)
  {
    if (list.empty())
    {
      o << "{}";
    }
    else
    {
      const bool is_hex = (o.flags() & std::ios_base::hex);
      o << "{ ";
      bool first = true;
      for (const auto& value : list)
      {
        if (!first)
        {
          o << ", ";
        }
        if (is_hex)
        {
          o << "0x";
        }
        o << value;
        first = false;
      }
      o << "\t- }";
    }
    return o;
  }

}  // namespace

TpcTimeFrameBuilderRun3::BcoMatchingInformation::BcoMatchingInformation(const std::string& name)
  : m_name(name)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_hNorm = new TH1D(TString(m_name.c_str()) + "_Normalization",  //
                     TString(m_name.c_str()) + " Normalization;Items;Count",
                     20, .5, 20.5);
  int i = 1;
  m_hNorm->GetXaxis()->SetBinLabel(i++, "SyncGTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "DC_STOP_SEND_GTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "HeartBeatGTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "HeartBeatFEE");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "HeartBeatFEEMatchedReference");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "HeartBeatFEEMatchedNew");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "HeartBeatFEEUnMatched");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "TriggerGTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "EnDATGTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "UnmatchedEnDATGTM");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FindGTMBCO");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FindGTMBCOMatchedExisting");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FindGTMBCOMatchedNew");
  m_hNorm->GetXaxis()->SetBinLabel(i++, "FindGTMBCOMatchedFailed");

  assert(i <= 20);
  m_hNorm->GetXaxis()->LabelsOption("v");
  hm->registerHisto(m_hNorm);

  m_hFEEClockAdjustment_MatchedReference = new TH1I(TString(m_name.c_str()) + "_FEEClockAdjustment_MatchedReference",  //
                                                    TString(m_name.c_str()) +
                                                        " FEEClockAdjustment for Matched Reference;Clock Adjustment [FEE Clock Cycle];Count",
                                                    512, -256 - .5, +256 - .5);
  hm->registerHisto(m_hFEEClockAdjustment_MatchedReference);
  m_hFEEClockAdjustment_MatchedNew = new TH1I(TString(m_name.c_str()) + "_FEEClockAdjustment_MatchedNew",  //
                                              TString(m_name.c_str()) +
                                                  " FEEClockAdjustment for Matched New;Clock Adjustment [FEE Clock Cycle];Count",
                                              512, -256 - .5, +256 - .5);
  hm->registerHisto(m_hFEEClockAdjustment_MatchedNew);

  m_hFEEClockAdjustment_Unmatched = new TH1I(TString(m_name.c_str()) + "_FEEClockAdjustment_Unmatched",  //
                                             TString(m_name.c_str()) +
                                                 " FEEClock Diff for unmatched;Clock Adjustment [FEE Clock Cycle];Count",
                                             512,
                                             -(1UL << m_FEE_CLOCK_BITS) - .5,
                                             +(1UL << m_FEE_CLOCK_BITS) - .5);
  hm->registerHisto(m_hFEEClockAdjustment_Unmatched);

  m_hGTMNewEventSpacing = new TH1I(TString(m_name.c_str()) +
                                       "_GTM_NewEventSpacing",  //
                                   TString(m_name.c_str()) +
                                       " Spacing between two events;Clock Diff [RHIC Clock Cycle];Count",
                                   1024, -.5, +1024 - .5);
  hm->registerHisto(m_hGTMNewEventSpacing);

  // m_hFindGTMBCO_MatchedExisting_BCODiff = new TH1I(TString(m_name.c_str()) + "_FindGTMBCO_MatchedExisting_BCODiff",  //
  //                                                  TString(m_name.c_str()) +
  //                                                      " find_gtm_bco matched to existing event clock diff;Clock Difference [FEE Clock Cycle];Count",
  //                                                  512, -256 - .5, +256 - .5);
  // hm->registerHisto(m_hFindGTMBCO_MatchedExisting_BCODiff);
  // m_hFindGTMBCO_MatchedNew_BCODiff = new TH1I(TString(m_name.c_str()) + "_FindGTMBCO_MatchedNew_BCODiff",  //
  //                                             TString(m_name.c_str()) +
  //                                                 " find_gtm_bco matched to new event clock diff;Clock Difference [FEE Clock Cycle];Count",
  //                                             512, -256 - .5, +256 - .5);
  // hm->registerHisto(m_hFindGTMBCO_MatchedNew_BCODiff);
}

//! whether reference bco has moved pass the given gtm_bco
bool TpcTimeFrameBuilderRun3::BcoMatchingInformation::isMoreDataRequired(const uint64_t& gtm_bco) const
{
  const uint64_t bco_correction = get_gtm_rollover_correction(gtm_bco);

  if (m_verbosity >= 2)
  {
    std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired entry"
              << " at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
              << " bco_correction = 0x" << std::hex << bco_correction << std::dec
              << std::endl;
  }

  if (m_bco_reference)
  {
    if (m_bco_reference.value().first > bco_correction + m_max_fee_sync_time)
    {
      if (m_verbosity >= 2)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                  << " at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                  << ". m_bco_reference.value().first = 0x" << std::hex << m_bco_reference.value().first << std::dec
                  << " bco_correction = 0x" << std::hex << bco_correction << std::dec
                  << ". satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                  << std::endl;
      }

      return false;
    }
  }

  if (!m_bco_heartbeat_list.empty())
  {
    if (m_bco_heartbeat_list.back().first > bco_correction + m_max_fee_sync_time)
    {
      if (m_verbosity >= 2)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                  << "at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                  << ". m_bco_heartbeat_list.back().first = 0x" << std::hex << m_bco_heartbeat_list.back().first << std::dec
                  << " bco_correction = 0x" << std::hex << bco_correction << std::dec
                  << ". satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                  << std::endl;
      }

      return false;
    }

    if (m_verbosity > 4)
    {
      std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                << "at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                << ". m_bco_heartbeat_list.back().first = 0x" << std::hex << m_bco_heartbeat_list.back().first << std::dec
                << " bco_correction = 0x" << std::hex << bco_correction << std::dec
                << ". not yet satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                << std::endl;
    }
  }

  if (m_verbosity > 3)
  {
    std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
              << "at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
              << " bco_correction = 0x" << std::hex << bco_correction << std::dec << ": more data required"
              << " as their is NO m_bco_reference nor m_bco_heartbeat_list"
              << std::endl;

    std::cout << "  m_gtm_bco_trigger_map:" << std::endl;
    for (const auto& trig : m_gtm_bco_trigger_map)
    {
      std::cout << " - 0x" << std::hex << trig.first << std::dec << "(Diff = " << trig.first - bco_correction << ") " << std::endl;
    }

    std::cout << "  m_bco_matching_list:" << std::endl;
    for (const auto& trig : m_bco_matching_list)
    {
      std::cout << " - 0x" << std::hex << trig.second << std::dec << "(Diff = " << trig.second - bco_correction << ") " << std::endl;
    }
  }
  return true;
}

//___________________________________________________
std::optional<uint32_t> TpcTimeFrameBuilderRun3::BcoMatchingInformation::get_predicted_fee_bco(uint64_t gtm_bco) const
{
  // check proper initialization
  if (!is_verified() || !m_bco_reference)
  {
    return std::nullopt;
  }

  // check whether it is within the same FEE clock rollover window based on the reference candidate list
  {
    uint64_t latest_reference_bco = (*m_bco_reference).first;
    if (! m_bco_heartbeat_list.empty())  
    {
        latest_reference_bco = m_bco_heartbeat_list.back().first;  // get the latest heartbeat bco
    }

    if (get_bco_diff(gtm_bco , latest_reference_bco)*m_clock_ratio_numerator 
    > ((1U << (m_FEE_CLOCK_BITS -1))) * m_clock_ratio_denominator)
    {
      if (m_verbosity >= 3)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::get_predicted_fee_bco -"
                  << " GTM bco 0x" << std::hex << gtm_bco << std::dec
                  << " is too far from the latest heartbeat bco 0x" << std::hex << latest_reference_bco << std::dec
                  << " get_bco_diff(gtm_bco , latest_reference_bco) =" << get_bco_diff(gtm_bco , latest_reference_bco)
                  << " > " << ((1U << (m_FEE_CLOCK_BITS -1))) * m_clock_ratio_denominator / m_clock_ratio_numerator
                  << ", cannot predict fee bco" << std::endl;
      }
      return std::nullopt;
    }
  }

  // get gtm bco difference with proper rollover accounting
  const auto& bco_reference = *m_bco_reference;
  const int64_t gtm_bco_difference = int64_t(gtm_bco) - int64_t(bco_reference.first);

  static_assert(m_clock_ratio_numerator > 0);
  static_assert(m_clock_ratio_denominator > 0);

  // convert to fee bco with the exact Run3 30/8 ratio, and truncate to 20 bits
  const int64_t fee_bco_predicted = int64_t(bco_reference.second) +
                                    (gtm_bco_difference * m_clock_ratio_numerator) / m_clock_ratio_denominator;
  return uint32_t(static_cast<uint64_t>(fee_bco_predicted) & 0xFFFFFU);
}

//___________________________________________________
void TpcTimeFrameBuilderRun3::BcoMatchingInformation::print_gtm_bco_information() const
{
  if (!m_gtm_bco_trig_list.empty())
  {
    std::cout
        << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::print_gtm_bco_information -"
        << "\t- m_gtm_bco_trig_list: " << std::hex << m_gtm_bco_trig_list << std::dec
        << std::endl;

    // also print predicted fee bco
    if (is_verified())
    {
      std::list<uint32_t> fee_bco_predicted_list;
      std::transform(
          m_gtm_bco_trig_list.begin(),
          m_gtm_bco_trig_list.end(),
          std::back_inserter(fee_bco_predicted_list),
          [this](const uint64_t& gtm_bco)
          { return get_predicted_fee_bco(gtm_bco).value(); });

      std::cout
          << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::print_gtm_bco_information -"
          << "\t- m_gtm_bco_trig_list fee predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
    }
  }

  std::cout << "\t m_gtm_bco_dc_read = " << std::hex
            << m_gtm_bco_dc_read.first << " -> 0x" << m_gtm_bco_dc_read.second
            << std::dec << std::endl;
}

uint64_t TpcTimeFrameBuilderRun3::BcoMatchingInformation::
    get_gtm_rollover_correction(const uint64_t& gtm_bco) const
{
  // start with 40bit clock, enforced
  uint64_t gtm_bco_corrected = gtm_bco & ((uint64_t(1) << m_GTM_CLOCK_BITS) - 1);

  if (!m_bco_reference)
  {
    return gtm_bco_corrected;
  }

  // get the last GTM clock roll over
  const uint64_t& last_bco = m_bco_reference.value().first;
  const uint64_t last_bco_rollover = last_bco &
                                     (std::numeric_limits<uint64_t>::max() << m_GTM_CLOCK_BITS);

  // use the roll over of the last GTM clock reading
  gtm_bco_corrected += last_bco_rollover;

  // check if the rollover has advanced
  if (gtm_bco_corrected + (uint64_t(1) << (m_GTM_CLOCK_BITS - 1)) < last_bco)
  {
    gtm_bco_corrected += uint64_t(1) << m_GTM_CLOCK_BITS;
  }

  return gtm_bco_corrected;
}

//___________________________________________________
void TpcTimeFrameBuilderRun3::BcoMatchingInformation::save_bx_counter_sync_observation(
    uint64_t bx_counter_sync_gtm_bco,
    uint64_t bco_reference_gtm_bco,
    const TpcTimeFrameBuilderRun3::BcoMatchingInformation::m_gtm_fee_bco_matching_pair_t& bco_reference)
{
  if (m_bx_counter_sync_observation_count >= kMaxBXCounterSyncObservations)
  {
    return;
  }

  BXCounterSyncObservation& observation = m_bx_counter_sync_observations[m_bx_counter_sync_observation_count];
  observation.bx_counter_sync_gtm_bco = bx_counter_sync_gtm_bco;
  observation.bco_reference_gtm_bco = bco_reference_gtm_bco;
  observation.m_bco_reference = bco_reference;
  ++m_bx_counter_sync_observation_count;
}

void TpcTimeFrameBuilderRun3::BcoMatchingInformation::save_gtm_bco_information(const TpcTimeFrameBuilderRun3::gtm_payload& gtm_tagger)
{
  // append gtm_bco from taggers in this event to packet-specific list of available lv1_bco

  // save level1 trigger bco
  const bool& is_lvl1 = gtm_tagger.is_lvl1;
  const bool& is_endat = gtm_tagger.is_endat;
  const bool& is_modebit = gtm_tagger.is_modebit;
  const uint64_t gtm_bco = get_gtm_rollover_correction(gtm_tagger.bco);

  if (is_lvl1)
  {
    assert(m_hNorm);
    m_hNorm->Fill("TriggerGTM", 1);

    assert(m_hGTMNewEventSpacing);
    if (!m_gtm_bco_trig_list.empty())
    {
      m_hGTMNewEventSpacing->Fill(gtm_bco - m_gtm_bco_trig_list.back());
    }
    m_gtm_bco_trig_list.push_back(gtm_bco);
  }

  // also save ENDDAT bco
  else if (is_endat)
  {
    assert(m_hNorm);
    m_hNorm->Fill("EnDATGTM", 1);

    // add to list if difference to last entry is big enough
    if (m_gtm_bco_trig_list.empty() || (gtm_bco - m_gtm_bco_trig_list.back()) > m_max_lv1_endat_bco_diff)
    {
      assert(m_hNorm);
      m_hNorm->Fill("UnmatchedEnDATGTM", 1);

      if (!m_gtm_bco_trig_list.empty())
      {
        assert(m_hGTMNewEventSpacing);
        m_hGTMNewEventSpacing->Fill(gtm_bco - m_gtm_bco_trig_list.back());
      }
      m_gtm_bco_trig_list.push_back(gtm_bco);
    }
  }

  // also save hearbeat bco
  else if (is_modebit)
  {
    // get modebits
    const uint64_t& modebits = gtm_tagger.modebits;
    if (modebits == ELINK_HEARTBEAT_T)
    {
      assert(m_hNorm);
      m_hNorm->Fill("HeartBeatGTM", 1);

      auto predicted_fee_bco = get_predicted_fee_bco(gtm_bco);
      if (predicted_fee_bco)
      {
        m_bco_heartbeat_list.emplace_back(gtm_bco, predicted_fee_bco.value());
      }
      else
      {
        if (m_verbosity > 1)
        {
          std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                    << "\t- Warning: predicted_fee_bco is not available for gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                    << ". Skipping heartbeat candidate." << std::endl;
        }
      }

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                  << "\t- found heartbeat candidate "
                  << "at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                  << ". Current m_bco_heartbeat_list:"
                  << std::endl;

        for (const m_gtm_fee_bco_matching_pair_t& bco : m_bco_heartbeat_list)
        {
          std::cout << "\t- gtm_bco = 0x" << std::hex << bco.first << std::dec
                    << "\t- fee_bco = 0x" << std::hex << bco.second << std::dec
                    << std::endl;
        }
      }

      while (m_bco_heartbeat_list.size() > m_max_bco_heartbeat_list_size)
      {
        if (m_verbosity > 1)
        {
          uint64_t bco = m_bco_heartbeat_list.begin()->first;
          std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_from_modebits"
                    << "Warning: m_bco_heartbeat_list is full"
                    << "\t- drop unprocessed heart beat in queue "
                    << "at gtm_bco = 0x" << std::hex << bco
                    << std::dec
                    << ". Unprocessed heartbeats in queue with size of " << m_bco_heartbeat_list.size()
                    << std::endl;
        }

        m_bco_heartbeat_list.pop_front();
      }

    }  //     if (modebits & (1U << ELINK_HEARTBEAT_T))

    if (modebits == BX_COUNTER_SYNC_T)  // initiate synchronization of clock sync
    {
      assert(m_hNorm);
      m_hNorm->Fill("SyncGTM", 1);

      // get BCO and assign
      const uint64_t bco_reference_gtm_bco = gtm_bco + kBXCounterSyncGtmBcoOffset;
      const m_gtm_fee_bco_matching_pair_t bx_counter_sync_reference =
          std::make_pair(bco_reference_gtm_bco, static_cast<uint32_t>(kBXCounterSyncFEEBcoOffset));
      m_verified_from_modebits = true;
      m_bco_reference = bx_counter_sync_reference;
      save_bx_counter_sync_observation(gtm_bco, bco_reference_gtm_bco, bx_counter_sync_reference);
      m_bco_heartbeat_list.clear();

      if (m_verbosity)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_from_modebits"
                  << "\t- found reference from modebits BX_COUNTER_SYNC_T "
                  << "at gtm_bco = 0x" << std::hex << gtm_bco
                  << " reference gtm_bco = 0x" << bco_reference_gtm_bco << std::dec
                  << " GTM sync offset = " << kBXCounterSyncGtmBcoOffset
                  << " FEE sync offset = " << kBXCounterSyncFEEBcoOffset
                  << std::endl;
      }
    }  //     if (modebits == BX_COUNTER_SYNC_T)  // initiate synchronization of clock sync

    if (modebits == DC_STOP_SEND_T)
    {
      assert(m_hNorm);
      m_hNorm->Fill("DC_STOP_SEND_GTM", 1);

      // save the gtm_bco for the digital current readout
      m_gtm_bco_dc_read.first = gtm_bco;
      if (is_verified())
      {
        m_gtm_bco_dc_read.second = get_predicted_fee_bco(gtm_bco).value();  // NOLINT(bugprone-unchecked-optional-access)
      }
      else
      {
        m_gtm_bco_dc_read.second = 0;  // not verified, so no reference clock sync available
      }

      if (m_verbosity > 2)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                  << "\t- found DC stop send modebit "
                  << "at gtm_bco = 0x" << std::hex << gtm_bco << std::dec
                  << std::endl;
      }
    }
  }
}

//___________________________________________________
std::optional<uint64_t> TpcTimeFrameBuilderRun3::BcoMatchingInformation::find_reference_heartbeat(const TpcTimeFrameBuilderRun3::fee_payload& HeartBeatPacket)
{
  assert(m_hNorm);
  m_hNorm->Fill("HeartBeatFEE", 1);

  // make sure the bco matching is properly initialized and historical valid
  if (!is_verified())
  {
    return std::nullopt;
  }

  assert(HeartBeatPacket.type == HEARTBEAT_T);
  const uint32_t& fee_bco = HeartBeatPacket.bx_timestamp;

  if (m_bco_reference)
  {
    const uint64_t& gtm_bco = m_bco_reference.value().first;
    const uint32_t& fee_bco_predicted = m_bco_reference.value().second;
    // check if the predicted fee bco matches the actual fee bco
    if (get_fee_bco_diff(fee_bco_predicted, fee_bco) < m_max_fee_bco_diff)
    {
      // Keep QA for matched heartbeat, but do not update clock reference from heartbeat.
      if (verbosity() > 1)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - found a matched reference heartbeat; clock reference update disabled: "
                  << std::hex
                  << "\t- fee_bco: 0x" << fee_bco
                  << "\t- predicted: 0x" << fee_bco_predicted
                  << "\t- gtm_bco: 0x" << gtm_bco
                  << std::dec
                  << std::endl;
      }

      assert(m_hFEEClockAdjustment_MatchedReference);
      m_hFEEClockAdjustment_MatchedReference->Fill(int64_t(fee_bco) - int64_t(fee_bco_predicted), 1);

      m_hNorm->Fill("HeartBeatFEEMatchedReference", 1);

      return gtm_bco;
    }
  }

  for (const m_gtm_fee_bco_matching_pair_t& bco : m_bco_heartbeat_list)
  {
    const uint64_t gtm_bco = bco.first;
    const uint32_t fee_bco_predicted = bco.second;

    // check if the predicted fee bco matches the actual fee bco
    if (get_fee_bco_diff(fee_bco_predicted, fee_bco) < m_max_fee_bco_diff)
    {
      if (verbosity() > 1)
      {
        std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - found a new reference candidate heartbeat; clock reference update disabled: "
                  << std::hex
                  << "\t- fee_bco: 0x" << fee_bco
                  << "\t- predicted: 0x" << fee_bco_predicted
                  << "\t- gtm_bco: 0x" << gtm_bco
                  << "\t- previous reference gtm_bco: 0x" << m_bco_reference.value().first   // NOLINT(bugprone-unchecked-optional-access)
                  << "\t- previous reference fee_bco: 0x" << m_bco_reference.value().second  // NOLINT(bugprone-unchecked-optional-access)
                  << std::dec
                  << std::endl;
      }
      // Keep QA for matched candidate heartbeat, but do not replace the clock reference or trim candidates.
      if (m_verbosity > 1)
      {
        std::cout << "\t- clock reference update from heartbeat is disabled; candidate list retained at size "
                  << m_bco_heartbeat_list.size() << std::endl;
      }

      assert(m_hFEEClockAdjustment_MatchedNew);
      m_hFEEClockAdjustment_MatchedNew->Fill(int64_t(fee_bco) - int64_t(fee_bco_predicted), 1);

      m_hNorm->Fill("HeartBeatFEEMatchedNew", 1);
      return gtm_bco;
    }

    if (verbosity() > 1)
    {
      std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - unmatched heartbeat: "
                << std::hex
                << "\t- fee_bco: 0x" << fee_bco
                << "\t- predicted: 0x" << fee_bco_predicted
                << "\t- gtm_bco: 0x" << gtm_bco
                << std::dec
                << std::endl;
    }

    assert(m_hFEEClockAdjustment_Unmatched);
    m_hFEEClockAdjustment_Unmatched->Fill(int64_t(fee_bco) - int64_t(fee_bco_predicted), 1);
  }  //   for (const auto& bco : m_bco_heartbeat_list)

  if (verbosity() > 1)
  {
    std::cout << "TpcTimeFrameBuilderRun3[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - WARNING: failed match for fee_bco = 0x" << std::hex << fee_bco << std::dec << std::endl;
  }
  m_hNorm->Fill("HeartBeatFEEUnMatched", 1);
  return std::nullopt;
}

//___________________________________________________
void TpcTimeFrameBuilderRun3::BcoMatchingInformation::cleanup()
{
  // remove old gtm_bco and matching
  while (m_gtm_bco_trig_list.size() > m_max_matching_data_size)
  {
    m_gtm_bco_trig_list.pop_front();
  }
  while (m_bco_matching_list.size() > m_max_matching_data_size)
  {
    m_bco_matching_list.pop_front();
  }
}

//___________________________________________________
void TpcTimeFrameBuilderRun3::BcoMatchingInformation::cleanup(uint64_t ref_bco)
{
  // erase all elements from bco_list that are less than or equal to ref_bco
  m_gtm_bco_trig_list.erase(std::remove_if(m_gtm_bco_trig_list.begin(), m_gtm_bco_trig_list.end(),
                                           [ref_bco](const uint64_t& bco)
                                           { return bco <= ref_bco; }),
                            m_gtm_bco_trig_list.end());

  // erase all elements from bco_list that are less than or equal to ref_bco
  m_bco_matching_list.erase(std::remove_if(m_bco_matching_list.begin(), m_bco_matching_list.end(),
                                           [ref_bco](const m_fee_gtm_bco_matching_pair_t& pair)
                                           {
                                             return pair.second <= ref_bco;
                                           }),
                            m_bco_matching_list.end());
}

void TpcTimeFrameBuilderRun3::fillBadFeeMap()
{
  const std::string filename = CDBInterface::instance()->getUrl("TPC_DECODER_BAD_FEE");

  if (filename.empty())
  {
    if (m_verbosity > 0)
    {
      std::cout << "TpcTimeFrameBuilderRun3::fillBadFeeMap - no file found for TPC_DECODER_BAD_FEE, not filling bad fee map" << std::endl;
    }
    return;
  }

  CDBTTree cdbtree(filename);
  cdbtree.LoadCalibrations();

  const int nentries = cdbtree.GetSingleIntValue("N_MASKED_FEES");

  for (int i = 0; i < nentries; i++)
  {
    m_maskedFEEs[cdbtree.GetIntValue(i, "EBDC")].insert(cdbtree.GetIntValue(i, "FEEID"));
  }
}
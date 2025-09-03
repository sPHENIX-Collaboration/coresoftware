#include "TpcTimeFrameBuilder.h"

#include <Event/oncsSubConstants.h>
#include <Event/packet.h>

#include <ffarawobjects/TpcRawHitv2.h>
#include <ffarawobjects/TpcRawHitv3.h>
#include <phool/PHTimer.h>  // for PHTimer

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <qautils/QAHistManagerDef.h>
#include <fun4all/PHTFileServer.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>

#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <tuple>   // For std::tie

using namespace std;

TpcTimeFrameBuilder::TpcTimeFrameBuilder(const int packet_id)
  : m_packet_id(packet_id)
  , m_HistoPrefix("TpcTimeFrameBuilder_Packet" + to_string(packet_id))
{
  for (int fee = 0; fee < MAX_FEECOUNT; ++fee)
  {
    m_bcoMatchingInformation_vec.emplace_back(
        std::string("BcoMatchingInformation_Packet") + to_string(packet_id) + "_FEE" + std::to_string(fee));
  }

  m_feeData.resize(MAX_FEECOUNT);

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_packetTimer = new PHTimer("TpcTimeFrameBuilder_Packet" + to_string(packet_id));

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  m_hNorm = new TH1D(TString(m_HistoPrefix.c_str()) + "_Normalization",  //
                     TString(m_HistoPrefix.c_str()) + " Normalization;Items;Count",
                     20, .5, 20.5);
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

  assert(i <= 20);
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

  h_ProcessPacket_Time = new TH2I(TString(m_HistoPrefix.c_str()) + "_ProcessPacket_Time",  //
                                  TString(m_HistoPrefix.c_str()) +
                                      " Time cost to run ProcessPacket();Call counts;Time elapsed per call [ms];Count",
                                  100, 0, 30e6, 100, 0, 10);
  hm->registerHisto(h_ProcessPacket_Time);
}

TpcTimeFrameBuilder::~TpcTimeFrameBuilder()
{
  for (auto& timeFrameEntry : m_timeFrameMap)
  {
    while (!timeFrameEntry.second.empty())
    {
      delete timeFrameEntry.second.back();
      timeFrameEntry.second.pop_back();
    }
  }

  if (m_packetTimer)
  {
    delete m_packetTimer;
  }

  if (m_digitalCurrentDebugTTree)
  {
    delete m_digitalCurrentDebugTTree;
  }
}

void TpcTimeFrameBuilder::setVerbosity(const int i)
{
  m_verbosity = i;

  for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    bcoMatchingInformation.set_verbosity(i);
  }
}

bool TpcTimeFrameBuilder::isMoreDataRequired(const uint64_t& gtm_bco) const
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

std::vector<TpcRawHit*>& TpcTimeFrameBuilder::getTimeFrame(const uint64_t& gtm_bco)
{
  assert(m_hNorm);
  uint64_t bclk_rollover_corrected = m_bcoMatchingInformation_vec[0].get_gtm_rollover_correction(gtm_bco);

  // cleanup old unused matching info after completion of a time frame
  for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    bcoMatchingInformation.cleanup(bclk_rollover_corrected);
  }

  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
              << ": getTimeFrame for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
              << ": bclk_rollover_corrected: 0x" << std::hex << bclk_rollover_corrected << std::dec
              << std::endl;
  }

  for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)
  {
    if (it->first + GL1_BCO_MATCH_WINDOW < bclk_rollover_corrected)
    {
      if (m_verbosity >= 2)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
                  << ":DROPPED: BCO " << std::hex << it->first << std::dec
                  << " dropped for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
                  << " and bclk_rollover_corrected 0x" << std::hex
                  << bclk_rollover_corrected << std::dec << ". m_timeFrameMap:" << std::endl;

        // if (m_verbosity >= 3)
        {
          for (const auto& timeframe : m_timeFrameMap)
          {
            std::cout << "- BCO in map: 0x" << std::hex << timeframe.first << std::dec
                      << "(Diff:" << int64_t(timeframe.first) - int64_t(bclk_rollover_corrected)
                      << ")"
                      << " size: " << timeframe.second.size()
                      << std::endl;
          }
        }
      }

      m_hNorm->Fill("GTM_TimeFrame_Dropped_Hit_Sum", it->second.size());
      assert(h_GTMClockDiff_Dropped);
      h_GTMClockDiff_Dropped->Fill(int64_t(it->first) - int64_t(bclk_rollover_corrected));
      for (const auto& hit : it->second)
      {
        delete hit;
      }
      it = m_timeFrameMap.erase(it);
    }
    else if (it->first < bclk_rollover_corrected + GL1_BCO_MATCH_WINDOW)
    {
      if (m_verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
                  << ":PASS: BCO " << std::hex << it->first << std::dec
                  << " matched for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
                  << " and bclk_rollover_corrected 0x" << std::hex
                  << bclk_rollover_corrected << std::dec << ". m_timeFrameMap:" << std::endl;

        for (const auto& timeframe : m_timeFrameMap)
        {
          std::cout << "- BCO in map: 0x" << std::hex << timeframe.first << std::dec
                    << "(Diff:" << int64_t(timeframe.first) - int64_t(bclk_rollover_corrected)
                    << ")"
                    << " size: " << timeframe.second.size()
                    << std::endl;
        }
      }

      m_hNorm->Fill("GTM_TimeFrame_Matched", 1);
      assert(h_GTMClockDiff_Matched);
      h_GTMClockDiff_Matched->Fill(int64_t(it->first) - int64_t(bclk_rollover_corrected));
      assert(h_TimeFrame_Matched_Size);
      h_TimeFrame_Matched_Size->Fill(it->second.size());
      m_hNorm->Fill("GTM_TimeFrame_Matched_Hit_Sum", it->second.size());
      m_UsedTimeFrameSet.push(it->first);
      return it->second;
    }
    else
    {
      break;
    }
  }

  if (m_verbosity >= 1)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
              << ":WARNING: BCO no match for gtm_bco: 0x" << std::hex << gtm_bco << std::dec
              << "and bclk_rollover_corrected 0x" << std::hex
              << bclk_rollover_corrected << std::dec << ". m_timeFrameMap:" << std::endl;

    if (m_verbosity >= 2)
    {
      for (const auto& timeframe : m_timeFrameMap)
      {
        std::cout << "- BCO in map: 0x" << std::hex << timeframe.first << std::dec
                  << "(Diff:" << int64_t(timeframe.first) - int64_t(bclk_rollover_corrected) << ")" << std::endl;
      }
    }
  }

  m_hNorm->Fill("GTM_TimeFrame_Unmatched", 1);
  assert(h_GTMClockDiff_Unmatched);
  for (const auto& timeframe : m_timeFrameMap)
  {
    h_GTMClockDiff_Unmatched->Fill(int64_t(timeframe.first) - int64_t(bclk_rollover_corrected));
  }
  static std::vector<TpcRawHit*> empty;
  return empty;
}

void TpcTimeFrameBuilder::CleanupUsedPackets(const uint64_t& bclk)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id << ": cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }

  while (! m_UsedTimeFrameSet.empty())
  {
    uint64_t bco_completed = m_UsedTimeFrameSet.front();
    m_UsedTimeFrameSet.pop();

    if (m_verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
                << ": cleaning up previous processed packet in m_timeFrameMap at clock 0x"
                << std::hex
                << bco_completed << std::dec
                << " for CleanupUsedPackets(const uint64_t& bclk) call at bclk 0x" << std::hex
                << bclk << std::dec
                << " Diff:" << int64_t(bco_completed) - int64_t(bclk) << std::endl;
    }

    auto it = m_timeFrameMap.find(bco_completed);

    // assert(it != m_timeFrameMap.end());  // strong workflow check disabled
    // this can happen if TPC GL1 tagger is shifted slight ahead of the GTM BCO, but within GL1_BCO_MATCH_WINDOW

    if (it != m_timeFrameMap.end())
    {
      while (!it->second.empty())
      {
        delete it->second.back();
        it->second.pop_back();
      }
      m_timeFrameMap.erase(it);
    }
  }

  uint64_t bclk_rollover_corrected = m_bcoMatchingInformation_vec[0].get_gtm_rollover_correction(bclk);

  assert(m_hFEEDataStream);

  for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)
  {
    if (it->first <= bclk_rollover_corrected)
    {
      int count = 0;
      while (!it->second.empty())
      {
        m_hFEEDataStream->Fill(it->second.back()->get_fee(), "HitUnusedBeforeCleanup", 1);
        delete it->second.back();
        it->second.pop_back();
        ++count;
      }

      if (m_verbosity >= 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "\t- packet " << m_packet_id
                  << ": cleaning up " << count << " TPC hits in m_timeFrameMap at clock 0x"
                  << std::hex
                  << it->first << std::dec
                  << " for <= bclk_rollover_corrected 0x" << std::hex
                  << bclk_rollover_corrected << std::dec
                  << " Diff:" << int64_t(it->first) - int64_t(bclk_rollover_corrected) << std::endl;
      }
      m_timeFrameMap.erase(it++);
    }
    else
    {
      break;
    }
  }  //   for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)
}

int TpcTimeFrameBuilder::ProcessPacket(Packet* packet)
{
  static size_t call_count = 0;
  ++call_count;

  if (m_verbosity > 1)
  {
    std::cout << "TpcTimeFrameBuilder::ProcessPacket: " << m_packet_id
              << "\t- Entry " << std::endl;
  }

  if (!packet)
  {
    cout << __PRETTY_FUNCTION__ << "\t- Error : Invalid packet, doing nothing" << endl;
    assert(packet);
    return 0;
  }

  if (m_hitFormat <0 )
  {
    if (packet->getHitFormat() != IDTPCFEEV4 && packet->getHitFormat() != IDTPCFEEV5 && packet->getHitFormat() != IDTPCFEEV6 )
    {
      cout << __PRETTY_FUNCTION__ << "\t- Error : expect packet format " << IDTPCFEEV4 
          << " or "<< IDTPCFEEV5<< " or "<< IDTPCFEEV6
          << "\t- but received packet format " << packet->getHitFormat() << ":" << endl;
      packet->identify();
      assert(packet->getHitFormat() == IDTPCFEEV4 || packet->getHitFormat() == IDTPCFEEV5 || packet->getHitFormat() == IDTPCFEEV6);
      return 0;
    }

    m_hitFormat = packet->getHitFormat();

    double clock_multiplier (0);
    if (m_hitFormat == IDTPCFEEV4)
    {
      // this is the clock multiplier from lvl1 to fee clock
      // Tested with Run24 data. Could be changable in future runs
      clock_multiplier = 4.262916255;
    }
    else if (m_hitFormat == IDTPCFEEV5 || packet->getHitFormat() == IDTPCFEEV6)
    {
      // version 46 FEE firmware 
      clock_multiplier = 30./8.;
    }
    else
    {
      assert(clock_multiplier>0);
    }

    if (m_verbosity >= 1)
    {
      cout << __PRETTY_FUNCTION__ << " set clock sync for hit format " << m_hitFormat 
      <<" with clock_multiplier = "<< clock_multiplier << endl;
    }
    for (BcoMatchingInformation& bcoMatchingInformation : m_bcoMatchingInformation_vec)
    {
      bcoMatchingInformation.set_gtm_clock_multiplier(clock_multiplier);
    }
  }
  else
  {
    if (packet->getHitFormat() != m_hitFormat)
    {
      cout << __PRETTY_FUNCTION__ << "\t- Error : expect the last packet format " << m_hitFormat
          << "\t- but received packet format " << packet->getHitFormat() << ":" << endl;
      packet->identify();
      assert((packet->getHitFormat() == m_hitFormat));
      return 0;
    }
  }
  assert((packet->getHitFormat() == m_hitFormat));

  if (m_packet_id != packet->getIdentifier())
  {
    cout << __PRETTY_FUNCTION__ << "\t- Error : mismatched packet with packet ID expectation of " << m_packet_id << ", but received";
    packet->identify();
    assert(m_packet_id == packet->getIdentifier());
    return 0;
  }

  assert(m_packetTimer);
  if ((m_verbosity == 1 && (call_count % 1000) == 0) || (m_verbosity > 1))
  {
    cout << __PRETTY_FUNCTION__ << "\t- : received packet ";
    packet->identify();

    m_packetTimer->print_stat();
  }
  m_packetTimer->restart();

  // //remove after testing
  // ;
  // cout <<"packet->lValue(0, N_TAGGER) = "<<packet->lValue(0, "N_TAGGER")<<endl;
  // cout <<"packet->lValue(0, NR_WF) = "<<packet->iValue(0, "NR_WF")<<endl;

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
    cout << __PRETTY_FUNCTION__ << "\t- : Warning : suspecious padding "
         << data_padding << "\t- in packet " << m_packet_id <<":" << endl;
    packet->identify();
    // packet->dump();
  }

  size_t dma_words_buffer = static_cast<size_t>(data_length) * 2 / DAM_DMA_WORD_LENGTH + 1;
  vector<dma_word> buffer(dma_words_buffer);

  int l2 = 0;
  packet->fillIntArray(reinterpret_cast<int*>(buffer.data()), data_length + DAM_DMA_WORD_LENGTH / 2, &l2, "DATA");
  
  if (data_padding != 0)
  {
    cout << __PRETTY_FUNCTION__ << "\t- :  data_length = " << data_length
         << "\t- data_padding = " << data_padding <<"\t l2 = "<<l2 << "\t- in packet " << m_packet_id <<":" << endl;
  }

  assert(l2 <= data_length);

  if(l2 < data_padding)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : Error : l2 from fillIntArray() is smaller than padding suggesting an invalid data: " << l2
         << "\t- in packet " << m_packet_id << ". Data length: " << data_length
         << ", data padding: " << data_padding <<". Ignore this packet: "<< endl;
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
    cout << __PRETTY_FUNCTION__ << "\t- : Warning : mismatch of RCDAQ data to DMA transfer. Dropping mismatched data: "
         << dma_residual << "\t- in packet " << m_packet_id << ". Dropping residual data : " << endl;

    assert(dma_words + 1 < buffer.size());
    const dma_word& last_dma_word_data = buffer[dma_words + 1];
    const uint16_t* last_dma_word = reinterpret_cast<const uint16_t*>(&last_dma_word_data);

    for (size_t i = 0; i < dma_residual; ++i)
    {
      cout << "\t- 0x" << hex << last_dma_word[i] << dec;
    }
    cout << endl;
  }

  if (m_verbosity > 1)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : packet" << m_packet_id << endl
         << "\t-   data_length = " << data_length << endl
         << "\t-   data_padding = " << data_padding << endl
         << "\t-   dma_words_buffer = " << dma_words_buffer << endl
         << "\t-   l2 = " << l2 << endl
         << "\t-   dma_words = " << dma_words << endl;
  }

  // demultiplexer
  for (size_t index = 0; index < dma_words; ++index)
  {
    const dma_word& dma_word_data = buffer[index];

    if (m_verbosity > 2)
    {
      cout << __PRETTY_FUNCTION__ << "\t- : processing DMA word "
           << index << "/" << dma_words << "\t- with header 0x"
           << hex << dma_word_data.dma_header << dec << endl;
    }

    if ((dma_word_data.dma_header & 0xFF00U) == FEE_MAGIC_KEY)
    {
      unsigned int fee_id = dma_word_data.dma_header & 0xffU;

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
        cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE ID " << fee_id << "\t- at position " << index << endl;
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
      cout << __PRETTY_FUNCTION__ << "\t- : Error : Unknown data type at position " << index << ": "
           << hex << buffer[index].dma_header << dec << endl;
      // not FEE data, e.g. GTM data or other stream, to be decoded
      m_hNorm->Fill("DMA_WORD_INVALID", 1);
    }
  }

  // sanity check for the timeframe size
  for (auto& timeframe : m_timeFrameMap)
  {
    if (timeframe.second.size() > kMaxRawHitLimit)
    {
      cout << __PRETTY_FUNCTION__ << "\t- : Warning : impossible amount of hits in the same timeframe at BCO "
           << timeframe.first << "\t- : " << timeframe.second.size() << ", limit is " << kMaxRawHitLimit
           << ". Dropping this time frame!"
           << endl;
      m_hNorm->Fill("TimeFrameSizeLimitError", 1);

      while (!timeframe.second.empty())
      {
        delete timeframe.second.back();
        timeframe.second.pop_back();
      }
    }
  }

  m_packetTimer->stop();
  assert(h_ProcessPacket_Time);
  h_ProcessPacket_Time->Fill(call_count, m_packetTimer->elapsed());

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTimeFrameBuilder::process_fee_data(unsigned int fee)
{
  assert(m_hFEEDataStream);

  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : processing FEE " << fee << "\t- with " << m_feeData[fee].size() << "\t- words" << endl;
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
        cout << __PRETTY_FUNCTION__ 
        << "\t- : processing FEE " << fee 
        << "\t- with digital packet" << endl;
      }
      
      m_hFEEDataStream->Fill(fee, "WordDigitalCurrentKeyWord", 1);
      is_digital_current = true;
    } //     if (data_buffer[3] == FEE_PACKET_MAGIC_KEY_3)
    else
    {

      if (data_buffer[1] != FEE_PACKET_MAGIC_KEY_1)
      {
        if (m_verbosity > 1)
        {
          cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE magic key at position 1 0x" << hex << data_buffer[1] << dec << endl;
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
          cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE magic key at position 2 0x" << hex << data_buffer[2] << dec << endl;
        }
        m_hFEEDataStream->Fill(fee, "WordSkipped", 1);
        data_buffer.pop_front();
        continue;
      }
      assert(data_buffer[2] == FEE_PACKET_MAGIC_KEY_2);

    }

    // valid packet
    const uint16_t & pkt_length = data_buffer[0];  // this is indeed the number of 10-bit words + 5 in this packet
    if (pkt_length > MAX_PACKET_LENGTH)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE pkt_length " << pkt_length << endl;
      }
      m_hFEEDataStream->Fill(fee, "InvalidLength", 1);
      data_buffer.pop_front();
      continue;
    }

    if (pkt_length + 1U > data_buffer.size())
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << "\t- : packet over buffer boundary for now, skip decoding and wait for more data: "
                                        " pkt_length = "
              << pkt_length
              << "\t- data_buffer.size() = " << data_buffer.size()
              << endl;
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

void TpcTimeFrameBuilder::process_fee_data_waveform(const unsigned int & fee, std::deque<uint16_t>& data_buffer)
{
  const uint16_t & pkt_length = data_buffer[0];

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

  if (! m_fastBCOSkip)
  {
    auto crc_parity = crc16_parity(fee, pkt_length);
    payload.calc_crc = crc_parity.first;
    payload.calc_parity = crc_parity.second;

    if (payload.data_crc != payload.calc_crc)
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << "\t- : CRC error in FEE "
              << fee << "\t- at position " << pkt_length - 1
              << ": data_crc = " << payload.data_crc
              << "\t- calc_crc = " << payload.calc_crc << endl;
      }
      m_hFEEDataStream->Fill(fee, "HitCRCError", 1);
      // continue;
    }

    if (payload.data_parity != payload.calc_parity)
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << "\t- : parity error in FEE "
              << fee << "\t- at position " << pkt_length - 1
              << ": data_parity = " << payload.data_parity
              << "\t- calc_parity = " << payload.calc_parity << endl;
      }
      m_hFEEDataStream->Fill(fee, "ParityError", 1);
      // continue;
    }
  }  //     if (not m_fastBCOSkip)

  assert(fee < m_bcoMatchingInformation_vec.size());
  BcoMatchingInformation& m_bcoMatchingInformation = m_bcoMatchingInformation_vec[fee];
  // gtm_bco matching
  if (payload.type == m_bcoMatchingInformation.HEARTBEAT_T)
  {
    if (m_verbosity > 1)
    {
      cout << __PRETTY_FUNCTION__
            << "\t- : received heartbeat packet from FEE " << fee << endl;
    }

    // if bco matching information is still not verified, drop the packet
    if (! m_bcoMatchingInformation.is_verified())
    {
      m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncUnavailable", 1);

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilder::process_fee_data - bco_matching not verified for heart beat, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
    }
    else  //       if (not m_bcoMatchingInformation.is_verified())
    {
      const optional<uint64_t> result = m_bcoMatchingInformation.find_reference_heartbeat(payload);
      m_hFEEDataStream->Fill(fee, "PacketHeartBeat", 1);

      if (result)
      {
        // assign gtm bco
        payload.gtm_bco = result.value();
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
  else if (! m_fastBCOSkip)  //     if (payload.type == m_bcoMatchingInformation.HEARTBEAT_T)
  {
    m_hFEEChannelPacketCount->Fill(fee * MAX_CHANNELS + payload.channel, 1);

    // if bco matching information is still not verified, drop the packet
    if (! m_bcoMatchingInformation.is_verified())
    {
      m_hFEEDataStream->Fill(fee, "PacketClockSyncUnavailable", 1);

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilder::process_fee_data - bco_matching not verified, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
    }
    else
    {
      const optional<uint64_t> result = m_bcoMatchingInformation.find_gtm_bco(payload.bx_timestamp);

      if (result)
      {
        // assign gtm bco
        payload.gtm_bco = result.value();
        m_hFEEDataStream->Fill(fee, "PacketClockSyncOK", 1);
      }
      else
      {
        if (m_verbosity > 1)
        {
          std::cout << "TpcTimeFrameBuilder::process_fee_data - WARNING: bco_matching failed!" << std::endl;
          m_bcoMatchingInformation.print_gtm_bco_information();
        }
        m_hFEEDataStream->Fill(fee, "PacketClockSyncError", 1);

        // skip the waverform
      }
    }
  }

  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : received data packet "
          << "\t- from FEE " << fee << endl
          << "\t- pkt_length = " << pkt_length << endl
          << "\t- type = " << payload.type << endl
          << "\t- adc_length = " << payload.adc_length << endl
          << "\t- sampa_address = " << payload.sampa_address << endl
          << "\t- sampa_channel = " << payload.sampa_channel << endl
          << "\t- channel = " << payload.channel << endl
          << "\t- bx_timestamp = 0x" << hex << payload.bx_timestamp << dec << endl
          << "\t- bco = 0x" << hex << payload.gtm_bco << dec << endl
          << "\t- data_crc = 0x" << hex << payload.data_crc << dec << endl
          << "\t- calc_crc = 0x" << hex << payload.calc_crc << dec << endl
          << "\t- data_parity = 0x" << hex << payload.data_parity << dec << endl
          << "\t- calc_parity = 0x" << hex << payload.calc_parity << dec << endl;
  }

  if ((! m_fastBCOSkip) && payload.gtm_bco > 0)
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
        cout << __PRETTY_FUNCTION__ << ": nsamp: " << nsamp
              << "+ pos: " << pos
              << " pkt_length: " << pkt_length << " start_t:" << start_t << endl;
      }

      if (pos + nsamp > pkt_length)
      {
        if (m_verbosity > 1)
        {
          cout << __PRETTY_FUNCTION__ << ": WARNING : nsamp: " << nsamp
                << "+ pos: " << pos
                << " > pkt_length: " << pkt_length << ", format error over length: " << endl;

          for (int print_pos = 0; print_pos <= pkt_length; ++print_pos)
          {
            cout << "\t[" << print_pos << "]=0x" << hex << data_buffer[print_pos] << dec << "(" << data_buffer[print_pos] << ")";
          }
          cout << endl;
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
        ++data_buffer_iterator;  //data_buffer[pos++];
      }
      payload.waveforms.emplace_back(start_t, std::move(adc));

      //   // an exception to deal with the last sample that is missing in the current hit format
      //   if (pos + 1 == pkt_length) break;
    }

    if (pos != pkt_length)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__ << ": WARNING : residual data at the end of decoding:"
              << " pos: " << pos
              << " <pkt_length: " << pkt_length << ", format error under length" << endl;
      }
      m_hFEEDataStream->Fill(fee, "HitFormatErrorMismatchedLength", 1);
    }

    // valid packet in the buffer, create a new hit
    if (payload.type != m_bcoMatchingInformation.HEARTBEAT_T)
    {
      TpcRawHitv3* hit = new TpcRawHitv3();
      m_timeFrameMap[payload.gtm_bco].push_back(hit);

      hit->set_bco(payload.bx_timestamp);
      hit->set_packetid(m_packet_id);
      hit->set_fee(fee);
      hit->set_channel(payload.channel);
      hit->set_type(payload.type);
      // hit->set_checksum(payload.data_crc);
      hit->set_checksumerror(payload.data_crc != payload.calc_crc);
      // hit->set_parity(payload.data_parity);
      hit->set_parityerror(payload.data_parity != payload.calc_parity);

      for (pair<uint16_t, std::vector<uint16_t>>& waveform : payload.waveforms)
      {
        hit->move_adc_waveform(waveform.first, std::move(waveform.second));
      }
    }
  }  //     if (not m_fastBCOSkip)

  return  ;
}

void TpcTimeFrameBuilder::process_fee_data_digital_current(const unsigned int & fee, std::deque<uint16_t>& data_buffer)
{
  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : processing digital_current data " << endl;
  }
  m_hFEEDataStream->Fill(fee, "DigitalCurrent", 1);
  const uint16_t & pkt_length = data_buffer[0];

  if (pkt_length != HEADER_LENGTH + digital_current_payload::MAX_CHANNELS * 2 * 2)
  {
    if (m_verbosity > 1)
    {
      cout << __PRETTY_FUNCTION__ << "\t- : Error : Invalid FEE pkt_length " << pkt_length
           << ", expected at least " << HEADER_LENGTH + digital_current_payload::MAX_CHANNELS * 2 * 2
           << endl;
    }
    m_hFEEDataStream->Fill(fee, "DigitalCurrentFormatErrorMismatchedLength", 1);
    return;
  }

  digital_current_payload payload;

  payload.fee           = fee;
  payload.pkt_length    = pkt_length;
  payload.sampa_address = (data_buffer[4] >> 5U) & 0xfU;
  // payload.sampa_max_channel = data_buffer[4] & 0x1fU;
  payload.channel       = data_buffer[4] & 0x1ffU;
  // payload.type          = data_buffer[3];
  payload.bx_timestamp  = ((data_buffer[6] & 0x3ffU) << 10U) | (data_buffer[5] & 0x3ff);

  uint16_t  pos = HEADER_LENGTH;
  for(int ich = 0; ich<digital_current_payload::MAX_CHANNELS; ich++)
  {
    payload.current[ich] = ((unsigned int)data_buffer[pos])<<16U | ((unsigned int)data_buffer[pos+1U]);
    pos++; pos++;
    payload.nsamples[ich] = ((unsigned int)data_buffer[pos])<<16U | ((unsigned int)data_buffer[pos+1U]);
    pos++; pos++;
  }

  if (pos != pkt_length)
  {
    if (m_verbosity> 1)
    {
      cout << __PRETTY_FUNCTION__ << "\t- : Warning : residual data at the end of decoding:"
           << " pos: " << pos
           << " <pkt_length: " << pkt_length << ", format error under length" << endl;
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
      cout << __PRETTY_FUNCTION__ << "\t- : CRC error in FEE "
           << fee << "\t- at position " << pkt_length - 1
           << ": data_crc = " << payload.data_crc
           << "\t- calc_crc = " << payload.calc_crc << endl;
    }
    m_hFEEDataStream->Fill(fee, "DigitalCurrentCRCError", 1);
    // continue;
  }

  assert(fee < m_bcoMatchingInformation_vec.size());
  BcoMatchingInformation& m_bcoMatchingInformation = m_bcoMatchingInformation_vec[fee];
  std::tie(payload.gtm_bco, payload.bx_timestamp_predicted) = m_bcoMatchingInformation.find_dc_read_bco();

  if (m_verbosity>2)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : received digital current packet "
         << "\t- from FEE " << fee << endl
         << "\t- pkt_length = " << pkt_length << endl
         << "\t- sampa_address = " << payload.sampa_address << endl
         << "\t- channel = " << payload.channel << endl
         << "\t- bx_timestamp = 0x" << hex << payload.bx_timestamp << dec << endl
         << "\t- gtm_bco = 0x" << hex << payload.gtm_bco << dec << endl
         << "\t- bx_timestamp_predicted = 0x" << hex << payload.bx_timestamp_predicted << dec << endl;

    cout << "\t- current:" ;
    for (int ich = 0; ich < digital_current_payload::MAX_CHANNELS; ich++)
    {
      cout << "\t[" << ich << "] = " << payload.current[ich]  ;
    }
    cout   << endl;
    cout << "\t- nsamples:" ;
    for (int ich = 0; ich < digital_current_payload::MAX_CHANNELS; ich++)
    {
      cout << "\t[" << ich << "] = " << payload.nsamples[ich] ;
    }
    cout   << endl;
    cout<< "\t- data_crc = 0x" << hex << payload.data_crc << dec << endl
       << "\t- calc_crc = 0x" << hex << payload.calc_crc << dec << endl;
  }

  if (m_digitalCurrentDebugTTree)
  {
    m_digitalCurrentDebugTTree->fill(payload);
  }

  return  ;
}

void TpcTimeFrameBuilder::SaveDigitalCurrentDebugTTree(const std::string &name)
{
  if (m_verbosity >= 1)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : Saving digital current debug TTree to " << name << endl;
  }

  m_digitalCurrentDebugTTree = new TpcTimeFrameBuilder::DigitalCurrentDebugTTree(name); 
}

TpcTimeFrameBuilder::DigitalCurrentDebugTTree::DigitalCurrentDebugTTree(const std::string &name)
: m_name(name) 
{  
  // open TFile
  PHTFileServer::get().open(m_name, "RECREATE");

  // cppcheck-suppress noCopyConstructor
  // cppcheck-suppress noOperatorEq
  m_tDigitalCurrent = new TTree("T_DigitalCurrent", "DigitalCurrent Debug TTree");
  assert(m_tDigitalCurrent);
  
  m_tDigitalCurrent->Branch("dc", &m_payload, 
    "gtm_bco/l:bx_timestamp_predicted/i:fee/s:pkt_length/s:channel/s:sampa_address/s:bx_timestamp/i:current[8]/i:nsamples[8]/i:data_crc/s:calc_crc/s");
}

TpcTimeFrameBuilder::DigitalCurrentDebugTTree::~DigitalCurrentDebugTTree()
{  
  // open TFile
  PHTFileServer::get().write(m_name);
}

void TpcTimeFrameBuilder::DigitalCurrentDebugTTree::fill(const TpcTimeFrameBuilder::digital_current_payload &payload)
{
  assert(m_tDigitalCurrent);

  m_payload = payload;
  m_tDigitalCurrent->Fill();
}

int TpcTimeFrameBuilder::decode_gtm_data(const TpcTimeFrameBuilder::dma_word& gtm_word)
{
  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << "\t- : processing GTM data " << endl;
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
    cout << __PRETTY_FUNCTION__ << "\t- GTM data : "
         << "\t- pkt_type = " << payload.pkt_type << endl
         << "\t- is_lvl1 = " << payload.is_lvl1 << endl
         << "\t- is_endat = " << payload.is_endat << endl
         << "\t- is_modebit = " << payload.is_modebit << endl
         << "\t- bco = 0x" << hex << payload.bco << dec << endl
         << "\t- lvl1_count = " << payload.lvl1_count << endl
         << "\t- endat_count = " << payload.endat_count << endl
         << "\t- last_bco = 0x" << hex << payload.last_bco << dec << endl
         << "\t- modebits =  0x" << hex << (int) payload.modebits << dec << endl
         << "\t- userbits =  0x" << hex << (int) payload.userbits << dec << endl;
  }

  if (payload.is_modebit)
  {
    if (payload.modebits == BcoMatchingInformation::ELINK_HEARTBEAT_T)
    {
      if (m_verbosity > 2)
      {
        cout << "\t- (Heartbeat modebit)" << endl;
      }
      assert(m_hNorm);
      m_hNorm->Fill("DMA_WORD_GTM_HEARTBEAT", 1);
    }

    if (payload.modebits == BcoMatchingInformation::DC_STOP_SEND_T)
    {
      if (m_verbosity > 2)
      {
        cout << "\t- (DC stop send modebit)" << endl;
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
        cout << __PRETTY_FUNCTION__ << "\t- : processing GTM data for FEE " << fee << endl;
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

uint16_t TpcTimeFrameBuilder::reverseBits(const uint16_t x) const
{
  uint16_t n = x;
  n = (static_cast<uint16_t>(n >> 1U) & 0x55555555U) | (static_cast<uint16_t>(n << 1U) & 0xaaaaaaaaU);
  n = (static_cast<uint16_t>(n >> 2U) & 0x33333333U) | (static_cast<uint16_t>(n << 2U) & 0xccccccccU);
  n = (static_cast<uint16_t>(n >> 4U) & 0x0f0f0f0fU) | (static_cast<uint16_t>(n << 4U) & 0xf0f0f0f0U);
  n = (static_cast<uint16_t>(n >> 8U) & 0x00ff00ffU) | (static_cast<uint16_t>(n << 8U) & 0xff00ff00U);
  // n = (n >> 16U) & 0x0000ffffU | (n << 16U) & 0xffff0000U;
  return n;
}

std::pair<uint16_t, uint16_t> TpcTimeFrameBuilder::crc16_parity(const uint32_t fee, const uint16_t l) const
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
  return make_pair(crc, data_parity);
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

TpcTimeFrameBuilder::BcoMatchingInformation::BcoMatchingInformation(const std::string& name)
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

  m_hFindGTMBCO_MatchedExisting_BCODiff = new TH1I(TString(m_name.c_str()) + "_FindGTMBCO_MatchedExisting_BCODiff",  //
                                                   TString(m_name.c_str()) +
                                                       " find_gtm_bco matched to existing event clock diff;Clock Difference [FEE Clock Cycle];Count",
                                                   512, -256 - .5, +256 - .5);
  hm->registerHisto(m_hFindGTMBCO_MatchedExisting_BCODiff);
  m_hFindGTMBCO_MatchedNew_BCODiff = new TH1I(TString(m_name.c_str()) + "_FindGTMBCO_MatchedNew_BCODiff",  //
                                              TString(m_name.c_str()) +
                                                  " find_gtm_bco matched to new event clock diff;Clock Difference [FEE Clock Cycle];Count",
                                              512, -256 - .5, +256 - .5);
  hm->registerHisto(m_hFindGTMBCO_MatchedNew_BCODiff);
}

//! whether reference bco has moved pass the given gtm_bco
bool TpcTimeFrameBuilder::BcoMatchingInformation::isMoreDataRequired(const uint64_t& gtm_bco) const
{
  const uint64_t bco_correction = get_gtm_rollover_correction(gtm_bco);

  if (m_verbosity>=2)
  {
    std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired entry"
              << " at gtm_bco = 0x" << hex << gtm_bco << dec
              << " bco_correction = 0x" << hex << bco_correction << dec
              << std::endl;
  }

  if (m_bco_reference)
  {
    if (m_bco_reference.value().first > bco_correction + m_max_fee_sync_time)
    {
      if (m_verbosity >= 2)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                  << " at gtm_bco = 0x" << hex << gtm_bco << dec
                  << ". m_bco_reference.value().first = 0x" << hex << m_bco_reference.value().first << dec
                  << " bco_correction = 0x" << hex << bco_correction << dec
                  << ". satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                  << std::endl;
      }

      return false;
    }
  }

  if (m_bco_reference_candidate_list.size() > 0)
  {
    if (m_bco_reference_candidate_list.back().first > bco_correction + m_max_fee_sync_time)
    {
      if (m_verbosity >= 2)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << ". m_bco_reference_candidate_list.back().first = 0x" << hex << m_bco_reference_candidate_list.back().first << dec
                  << " bco_correction = 0x" << hex << bco_correction << dec
                  << ". satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                  << std::endl;
      }

      return false;
    }
    else
    {
      if (m_verbosity > 4)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << ". m_bco_reference_candidate_list.back().first = 0x" << hex << m_bco_reference_candidate_list.back().first << dec
                  << " bco_correction = 0x" << hex << bco_correction << dec
                  << ". not yet satisified m_max_fee_sync_time = " << m_max_fee_sync_time
                  << std::endl;
      }
    }
  }

  if (m_verbosity > 3)
  {
    std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::isMoreDataRequired"
              << "at gtm_bco = 0x" << hex << gtm_bco << dec
              << " bco_correction = 0x" << hex << bco_correction << dec << ": more data required"
              << " as their is NO m_bco_reference nor m_bco_reference_candidate_list"
              << std::endl;

    std::cout << "  m_gtm_bco_trigger_map:" << std::endl;
    for (const auto& trig : m_gtm_bco_trigger_map)
    {
      std::cout << " - 0x" << hex << trig.first << dec << "(Diff = " << trig.first - bco_correction << ") " << std::endl;
    }

    std::cout << "  m_bco_matching_list:" << std::endl;
    for (const auto& trig : m_bco_matching_list)
    {
      std::cout << " - 0x" << hex << trig.second << dec << "(Diff = " << trig.second - bco_correction << ") " << std::endl;
    }
  }
  return true;
}

//___________________________________________________
std::optional<uint32_t> TpcTimeFrameBuilder::BcoMatchingInformation::get_predicted_fee_bco(uint64_t gtm_bco) const
{
  // check proper initialization
  if (!is_verified())
  {
    return std::nullopt;
  }

  // get gtm bco difference with proper rollover accounting
  const int64_t gtm_bco_difference = int64_t(gtm_bco) - int64_t(m_bco_reference.value().first);

  assert(m_multiplier>0);

  // convert to fee bco, and truncate to 20 bits
  const int64_t fee_bco_predicted = int64_t(m_bco_reference.value().second) + int64_t(m_multiplier * gtm_bco_difference);
  return uint32_t(static_cast<uint64_t>(fee_bco_predicted) & 0xFFFFFU);
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information() const
{
  if (!m_gtm_bco_trig_list.empty())
  {
    std::cout
        << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::print_gtm_bco_information -"
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
          [this](const uint64_t& gtm_bco) { return get_predicted_fee_bco(gtm_bco).value(); });

      std::cout
          << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::print_gtm_bco_information -"
          << "\t- m_gtm_bco_trig_list fee predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
    }
  }

  std::cout <<"\t m_gtm_bco_dc_read = " << std::hex 
    << m_gtm_bco_dc_read.first <<" -> 0x" << m_gtm_bco_dc_read.second 
    << std::dec << std::endl;
}

uint64_t TpcTimeFrameBuilder::BcoMatchingInformation::
    get_gtm_rollover_correction(const uint64_t& gtm_bco) const
{
  // start with 40bit clock, enforced
  uint64_t gtm_bco_corrected = gtm_bco & ((uint64_t(1) << m_GTM_CLOCK_BITS) - 1);

  if (! m_bco_reference)
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
void TpcTimeFrameBuilder::BcoMatchingInformation::save_gtm_bco_information(const TpcTimeFrameBuilder::gtm_payload& gtm_tagger)
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
    if (! m_gtm_bco_trig_list.empty())
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

      if (! m_gtm_bco_trig_list.empty())
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
        m_bco_reference_candidate_list.emplace_back(gtm_bco, predicted_fee_bco.value());
      }
      else
      {
        if (m_verbosity > 1)
        {
          std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                    << "\t- Warning: predicted_fee_bco is not available for gtm_bco = 0x" << hex << gtm_bco << dec
                    << ". Skipping heartbeat candidate." << std::endl;
        }
      }

      if (m_verbosity > 1)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                  << "\t- found heartbeat candidate "
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << ". Current m_bco_reference_candidate_list:"
                  << std::endl;

        for (const m_gtm_fee_bco_matching_pair_t& bco : m_bco_reference_candidate_list)
        {
          std::cout << "\t- gtm_bco = 0x" << hex << bco.first << dec
                    << "\t- fee_bco = 0x" << hex << bco.second << dec
                    << std::endl;
        }
      }

      while (m_bco_reference_candidate_list.size() > m_max_bco_reference_candidate_list_size)
      {
        if (m_verbosity > 1)
        {
          uint64_t bco = m_bco_reference_candidate_list.begin()->first;
          std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_from_modebits"
                    << "Warning: m_bco_reference_candidate_list is full"
                    << "\t- drop unprocessed heart beat in queue "
                    << "at gtm_bco = 0x" << hex << bco
                    << dec
                    << ". Unprocessed heartbeats in queue with size of " << m_bco_reference_candidate_list.size()
                    << std::endl;
        }

        m_bco_reference_candidate_list.pop_front();
      }

    }  //     if (modebits & (1U << ELINK_HEARTBEAT_T))

    if (modebits == BX_COUNTER_SYNC_T)  // initiate synchronization of clock sync
    {
      assert(m_hNorm);
      m_hNorm->Fill("SyncGTM", 1);

      // get BCO and assign
      m_verified_from_modebits = true;
      m_bco_reference = make_pair(gtm_bco, 0);
      m_bco_reference_candidate_list.clear();

      if (m_verbosity)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_from_modebits"
                  << "\t- found reference from modebits BX_COUNTER_SYNC_T "
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << std::endl;
      }
    } //     if (modebits == BX_COUNTER_SYNC_T)  // initiate synchronization of clock sync

    if (modebits == DC_STOP_SEND_T)
    {
      assert(m_hNorm);
      m_hNorm->Fill("DC_STOP_SEND_GTM", 1);

      // save the gtm_bco for the digital current readout
      m_gtm_bco_dc_read.first = gtm_bco;
      if (is_verified())
      {
        m_gtm_bco_dc_read.second = get_predicted_fee_bco(gtm_bco).value();
      }
      else
      {
        m_gtm_bco_dc_read.second = 0;  // not verified, so no reference clock sync available
      }

      if (m_verbosity > 2)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::save_gtm_bco_information"
                  << "\t- found DC stop send modebit "
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << std::endl;
      }
    }
  }
}

//___________________________________________________
std::optional<uint64_t> TpcTimeFrameBuilder::BcoMatchingInformation::find_reference_heartbeat(const TpcTimeFrameBuilder::fee_payload& HeartBeatPacket)
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
      // assign gtm bco
      m_bco_reference.value().second = fee_bco;

      if (verbosity() > 1)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - found an updated reference heartbeat and updated reference clock sync: "
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

  for (const m_gtm_fee_bco_matching_pair_t& bco : m_bco_reference_candidate_list)
  {
    const uint64_t gtm_bco = bco.first;
    const uint32_t fee_bco_predicted = bco.second;

    // check if the predicted fee bco matches the actual fee bco
    if (get_fee_bco_diff(fee_bco_predicted, fee_bco) < m_max_fee_bco_diff)
    {
      if (verbosity() > 1)
      {
        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - found a new reference canidate heartbeat and replaced reference clock sync: "
                  << std::hex
                  << "\t- fee_bco: 0x" << fee_bco
                  << "\t- predicted: 0x" << fee_bco_predicted
                  << "\t- gtm_bco: 0x" << gtm_bco
                  << "\t- previous reference gtm_bco: 0x" << m_bco_reference.value().first
                  << "\t- previous reference fee_bco: 0x" << m_bco_reference.value().second
                  << std::dec
                  << std::endl;
      }
      // assign gtm bco
      m_bco_reference = make_pair(gtm_bco, fee_bco);

      if (m_verbosity > 1)
      {
        std::cout << "\t- trimming m_bco_reference_candidate_list from size " << m_bco_reference_candidate_list.size() << std::endl;

        for (const m_gtm_fee_bco_matching_pair_t& bco_tmp : m_bco_reference_candidate_list)
        {
          std::cout << "\t\t- gtm_bco = 0x" << hex << bco_tmp.first << dec
                    << "\t\t- fee_bco = 0x" << hex << bco_tmp.second << dec
                    << std::endl;
        }
      }

      // remove the older candidate from the list
      while (m_bco_reference_candidate_list.begin()->first != gtm_bco)
      {
        m_bco_reference_candidate_list.pop_front();
      }
      m_bco_reference_candidate_list.pop_front();

      if (m_verbosity > 1)
      {
        std::cout << "\t- to size " << m_bco_reference_candidate_list.size() << std::endl;

        for (const m_gtm_fee_bco_matching_pair_t& bco_tmp : m_bco_reference_candidate_list)
        {
          std::cout << "\t\t- gtm_bco = 0x" << hex << bco_tmp.first << dec
                    << "\t\t- fee_bco = 0x" << hex << bco_tmp.second << dec
                    << std::endl;
        }
      }

      assert(m_hFEEClockAdjustment_MatchedNew);
      m_hFEEClockAdjustment_MatchedNew->Fill(int64_t(fee_bco) - int64_t(fee_bco_predicted), 1);

      m_hNorm->Fill("HeartBeatFEEMatchedNew", 1);
      return gtm_bco;
    }

    if (verbosity() > 1)
    {
      std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - unmatched heartbeat: "
                << std::hex
                << "\t- fee_bco: 0x" << fee_bco
                << "\t- predicted: 0x" << fee_bco_predicted
                << "\t- gtm_bco: 0x" << gtm_bco
                << std::dec
                << std::endl;
    }

    assert(m_hFEEClockAdjustment_Unmatched);
    m_hFEEClockAdjustment_Unmatched->Fill(int64_t(fee_bco) - int64_t(fee_bco_predicted), 1);
  }  //   for (const auto& bco : m_bco_reference_candidate_list)

  if (verbosity() > 1)
  {
    std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_reference_heartbeat - WARNING: failed match for fee_bco = 0x" << hex << fee_bco << dec << std::endl;
  }
  m_hNorm->Fill("HeartBeatFEEUnMatched", 1);
  return std::nullopt;
}

//___________________________________________________
std::optional<uint64_t> TpcTimeFrameBuilder::BcoMatchingInformation::find_gtm_bco(uint32_t fee_bco)
{
  if (verbosity() > 5)
  {
    std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_gtm_bco - entry: "
              << std::hex
              << "\t- fee_bco: 0x" << fee_bco
              << std::dec
              << "\t- is_verified(): " << (is_verified() ? "true" : "false")
              << std::endl;
  }

  // make sure the bco matching is properly initialized
  if (!is_verified())
  {
    return std::nullopt;
  }

  assert(m_hNorm);
  m_hNorm->Fill("FindGTMBCO", 1);

  // find matching gtm bco in map
  const auto bco_matching_iter = std::find_if(
      m_bco_matching_list.begin(),
      m_bco_matching_list.end(),
      [fee_bco](const m_fee_gtm_bco_matching_pair_t& pair) { return get_fee_bco_diff(pair.first, fee_bco) < m_max_fee_bco_diff; });

  if (bco_matching_iter != m_bco_matching_list.end())
  {
    m_hNorm->Fill("FindGTMBCOMatchedExisting", 1);
    assert(m_hFindGTMBCO_MatchedExisting_BCODiff);
    m_hFindGTMBCO_MatchedExisting_BCODiff->Fill(int64_t(fee_bco) - int64_t(bco_matching_iter->first));

    if (verbosity() > 3)
    {
      std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_gtm_bco - found existing FEE BCO: "
                << std::hex
                << "\t- fee_bco: 0x" << fee_bco
                << "\t- predicted: 0x" << bco_matching_iter->first
                << "\t- gtm_bco: 0x" << bco_matching_iter->second
                << std::dec
                << std::endl;
    }

    return bco_matching_iter->second;
  }
  else
  {
    // find element for which predicted fee_bco matches fee_bco, within limit
    const auto iter = std::find_if(
        m_gtm_bco_trig_list.begin(),
        m_gtm_bco_trig_list.end(),
        [this, fee_bco](const uint64_t& gtm_bco) { return get_fee_bco_diff(get_predicted_fee_bco(gtm_bco).value(), fee_bco) < m_max_gtm_bco_diff; });

    // check
    if (iter != m_gtm_bco_trig_list.end())
    {
      const uint64_t gtm_bco = *iter;

      m_hNorm->Fill("FindGTMBCOMatchedNew", 1);
      assert(m_hFindGTMBCO_MatchedNew_BCODiff);
      m_hFindGTMBCO_MatchedNew_BCODiff->Fill(int64_t(fee_bco) - int64_t(gtm_bco));

      if (verbosity() > 2)
      {
        const uint32_t fee_bco_predicted = get_predicted_fee_bco(gtm_bco).value();
        const uint32_t fee_bco_diff = get_bco_diff(fee_bco_predicted, fee_bco);

        std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_gtm_bco - new GL1 match: "
                  << std::hex
                  << "\t- fee_bco: 0x" << fee_bco
                  << "\t- predicted: 0x" << fee_bco_predicted
                  << "\t- gtm_bco: 0x" << gtm_bco
                  << std::dec
                  << "\t- difference: " << fee_bco_diff
                  << std::endl;
      }

      // save fee_bco and gtm_bco matching in map
      m_bco_matching_list.emplace_back(fee_bco, gtm_bco);

      // remove gtm bco from runing list
      m_gtm_bco_trig_list.erase(iter);

      // // update clock adjustment not applied for non HEARTBEAT_T
      // update_multiplier_adjustment(gtm_bco, fee_bco);

      return gtm_bco;
    }
    else
    {
      m_hNorm->Fill("FindGTMBCOMatchedFailed", 1);

      bool new_orphan = m_orphans.insert(fee_bco).second;

      if ((new_orphan && verbosity()) || (verbosity() > 3))
      {
        // find element for which predicted fee_bco is the closest to request
        const auto iter2 = std::min_element(
            m_gtm_bco_trig_list.begin(),
            m_gtm_bco_trig_list.end(),
            [this, fee_bco](const uint64_t& first, const uint64_t& second) { return get_bco_diff(get_predicted_fee_bco(first).value(), fee_bco) < get_bco_diff(get_predicted_fee_bco(second).value(), fee_bco); });

        const int fee_bco_diff = (iter2 != m_gtm_bco_trig_list.end()) ? get_bco_diff(get_predicted_fee_bco(*iter2).value(), fee_bco) : -1;

        if (m_verbosity >=2)
        {
          std::cout << "TpcTimeFrameBuilder[" << m_name << "]::BcoMatchingInformation::find_gtm_bco - match failed!"
                    << std::hex
                    << "\t- fee_bco: 0x" << fee_bco
                    << std::dec
                    << "\t- gtm_bco: 0x" << *iter2
                    << "\t- difference: " << fee_bco_diff
                    << std::endl;
        }
      }  //       if ((new_orphan and verbosity()) or (verbosity()>3))

      if (verbosity() > 3)
      {
        std::cout << "\t- m_gtm_bco_trig_list : " << std::endl;
        for (const auto& gtm_bco : m_gtm_bco_trig_list)
        {
          std::cout << "\t\t- 0x" << hex << gtm_bco << " -> 0x" << get_predicted_fee_bco(gtm_bco).value() << dec << std::endl;
        }

        std::cout << "\t- m_bco_matching_list : " << std::endl;
        for (const auto& iter_m_bco_matching_list : m_bco_matching_list)
        {
          std::cout << "\t\t- 0x" << hex << iter_m_bco_matching_list.first << " -> 0x" << iter_m_bco_matching_list.second << dec << std::endl;
        }

      }  //       if (verbosity()>3)

      return std::nullopt;
    }  // else
  }

  // never reached
  return std::nullopt;
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::cleanup()
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

  // clear orphans
  m_orphans.clear();
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::cleanup(uint64_t ref_bco)
{
  // erase all elements from bco_list that are less than or equal to ref_bco
  m_gtm_bco_trig_list.erase(std::remove_if(m_gtm_bco_trig_list.begin(), m_gtm_bco_trig_list.end(),
                                           [ref_bco](const uint64_t& bco) { return bco <= ref_bco; }),
                            m_gtm_bco_trig_list.end());

  // erase all elements from bco_list that are less than or equal to ref_bco
  m_bco_matching_list.erase(std::remove_if(m_bco_matching_list.begin(), m_bco_matching_list.end(),
                                           [ref_bco](const m_fee_gtm_bco_matching_pair_t& pair) {
                                             return pair.second <= ref_bco;
                                           }),
                            m_bco_matching_list.end());

  // clear orphans
  m_orphans.clear();
}

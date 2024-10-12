#include "TpcTimeFrameBuilder.h"

#include <Event/oncsSubConstants.h>
#include <Event/packet.h>

#include <ffarawobjects/TpcRawHitv2.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <qautils/QAHistManagerDef.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <stdint.h>
#include <cassert>
#include <memory>
#include <string>

using namespace std;

TpcTimeFrameBuilder::TpcTimeFrameBuilder(const int packet_id)
  : m_packet_id(packet_id)
  , m_HistoPrefix("TpcTimeFrameBuilder_Packet" + to_string(packet_id))
  , m_bcoMatchingInformation_vec(MAX_FEECOUNT)
{
  m_feeData.resize(MAX_FEECOUNT);

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1* h = new TH1D(TString(m_HistoPrefix.c_str()) + "_Normalization",  //
                    TString(m_HistoPrefix.c_str()) + " Normalization;Items;Count",
                    20, .5, 20.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Packet");
  h->GetXaxis()->SetBinLabel(i++, "Lv1-Taggers");
  h->GetXaxis()->SetBinLabel(i++, "EnDat-Taggers");
  h->GetXaxis()->SetBinLabel(i++, "ChannelPackets");
  h->GetXaxis()->SetBinLabel(i++, "Waveforms");

  h->GetXaxis()->SetBinLabel(i++, "DMA_WORD_GTM");
  h->GetXaxis()->SetBinLabel(i++, "DMA_WORD_FEE");
  h->GetXaxis()->SetBinLabel(i++, "DMA_WORD_FEE_INVALID");
  h->GetXaxis()->SetBinLabel(i++, "DMA_WORD_INVALID");

  h->GetXaxis()->SetBinLabel(i++, "TimeFrameSizeLimitError");

  assert(i <= 20);
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

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
                              MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5, 20, .5, 20.5);
  i = 1;
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordValid");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordSkipped");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "InvalidLength");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "RawHit");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitFormatError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "MissingLastHit");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitCRCError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitUnusedBeforeCleanup");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeat");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncUnavailable");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketHeartBeatClockSyncOK");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncUnavailable");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "PacketClockSyncOK");
  assert(i <= 20);
  hm->registerHisto(m_hFEEDataStream);



  m_hFEEClockAdjustment = new TH2I(TString(m_HistoPrefix.c_str()) + "_FEE_ClockAdjustment",  //
                              TString(m_HistoPrefix.c_str()) +
                                  " FEE Clock Adjustment;FEE ID;Type;Count",
                              MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5, 400, -1e-4, 1e-4);
  hm->registerHisto(m_hFEEClockAdjustment);

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
}

TpcTimeFrameBuilder::~TpcTimeFrameBuilder()
{
  for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end(); ++it)
  {
    while (!it->second.empty())
    {
      delete it->second.back();
      it->second.pop_back();
    }
  }
}

void TpcTimeFrameBuilder::setVerbosity(const int i)
{
  m_verbosity = i;

  for (auto& bcoMatchingInformation : m_bcoMatchingInformation_vec)
    bcoMatchingInformation.set_verbosity(i);
}

bool TpcTimeFrameBuilder::isMoreDataRequired() const
{
  // if (m_gtmData.size() == 0) return true;

  // if (m_gtmData.rbegin()->first - m_gtmData.begin()->first > kFEEDataTransmissionWindow)
  //   return false;

  return true;
}

void TpcTimeFrameBuilder::CleanupUsedPackets(const uint64_t bclk)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << " packet " << m_packet_id << ": cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }

  // uint64_t bclk_rollover_corrected = bclk;
  // if (m_gtmData.begin() != m_gtmData.end())
  // {
  //   if ((m_gtmData.begin()->first & GTMBCOmask_ValidBits) > bclk + (1ULL << (GTMBCObits - 1)))
  //   {
  //     bclk_rollover_corrected = (((m_gtmData.begin()->first >> GTMBCObits) + 1) << GTMBCObits) | (bclk & GTMBCOmask_ValidBits);
  //   }
  //   else if ((m_gtmData.begin()->first & GTMBCOmask_ValidBits) + (1ULL << (GTMBCObits - 1)) < bclk)
  //   {
  //     bclk_rollover_corrected = (((m_gtmData.begin()->first >> GTMBCObits) - 1) << GTMBCObits) | (bclk & GTMBCOmask_ValidBits);
  //   }
  //   else
  //   {
  //     bclk_rollover_corrected = (m_gtmData.begin()->first & GTMBCOmask_RollOverCounts) | (bclk & GTMBCOmask_ValidBits);
  //   }
  // }

  // assert(m_hFEEDataStream);

  // for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)
  // {
  //   if (it->first <= bclk_rollover_corrected)
  //   {
  //     while (!it->second.empty())
  //     {
  //       m_hFEEDataStream->Fill(it->second.back()->get_fee(), "HitUnusedBeforeCleanup", 1);
  //       delete it->second.back();
  //       it->second.pop_back();
  //     }
  //     m_timeFrameMap.erase(it++);
  //   }
  //   else
  //   {
  //     break;
  //   }
  // }  //   for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)

  // for (auto it = m_gtmData.begin(); it != m_gtmData.end();)
  // {
  //   if (it->first <= bclk)
  //   {
  //     m_gtmData.erase(it++);
  //   }
  //   else
  //   {
  //     break;
  //   }
  // }  //   for (auto it = m_gtmData.begin(); it != m_gtmData.end();)
}

int TpcTimeFrameBuilder::ProcessPacket(Packet* packet)
{
  if (m_verbosity > 1)
  {
    std::cout << "TpcTimeFrameBuilder::ProcessPacket: " << m_packet_id
              << " Entry " << std::endl;
  }

  if (!packet)
  {
    cout << __PRETTY_FUNCTION__ << " Error : Invalid packet, doing nothing" << endl;
    assert(packet);
    return 0;
  }

  if (packet->getHitFormat() != IDTPCFEEV4)
  {
    cout << __PRETTY_FUNCTION__ << " Error : expect packet format " << IDTPCFEEV4
         << " but received packet format " << packet->getHitFormat() << ":" << endl;
    packet->identify();
    assert(packet->getHitFormat() == IDTPCFEEV4);
    return 0;
  }

  if (m_packet_id != packet->getIdentifier())
  {
    cout << __PRETTY_FUNCTION__ << " Error : mismatched packet with packet ID expectation of " << m_packet_id << ", but received";
    packet->identify();
    assert(m_packet_id == packet->getIdentifier());
    return 0;
  }

  if (m_verbosity > 1)
  {
    cout << __PRETTY_FUNCTION__ << " : received packet ";
    packet->identify();
  }

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto(
      m_HistoPrefix + "_Normalization"));
  assert(h_norm);
  h_norm->Fill("Packet", 1);

  int data_length = packet->getDataLength();  // 32bit length
  assert(h_PacketLength);
  h_PacketLength->Fill(data_length);

  int data_padding = packet->getPadding();  // 32bit padding
  assert(h_PacketLength_Padding);
  h_PacketLength_Padding->Fill(data_padding);
  if (data_padding != 0)
  {
    cout << __PRETTY_FUNCTION__ << " : Warning : suspecious padding "
         << data_padding << " in packet " << m_packet_id << endl;
  }

  size_t dma_words_buffer = data_length * 2 / DAM_DMA_WORD_LENGTH + 1;
  vector<dma_word> buffer(dma_words_buffer);

  int l2 = 0;
  packet->fillIntArray(reinterpret_cast<int*>(buffer.data()), data_length + DAM_DMA_WORD_LENGTH / 2, &l2, "DATA");
  assert(l2 <= data_length);
  l2 -= data_padding;
  assert(l2 >= 0);

  size_t dma_words = l2 * 2 / DAM_DMA_WORD_LENGTH;
  size_t dma_residual = (l2 * 2) % DAM_DMA_WORD_LENGTH;
  assert(dma_words <= buffer.size());
  assert(h_PacketLength_Residual);
  h_PacketLength_Residual->Fill(dma_residual);
  if (dma_residual > 0)
  {
    cout << __PRETTY_FUNCTION__ << " : Warning : mismatch of RCDAQ data to DMA transfer. Dropping mismatched data: "
         << dma_residual << " in packet " << m_packet_id << ". Dropping residual data : " << endl;

    assert(dma_words + 1 < buffer.size());
    const dma_word& last_dma_word_data = buffer[dma_words + 1];
    const uint16_t* last_dma_word = reinterpret_cast<const uint16_t*>(&last_dma_word_data);

    for (size_t i = 0; i < dma_residual; ++i)
    {
      cout << " 0x" << hex << last_dma_word[i] << dec;
    }
    cout << endl;
  }

  if (m_verbosity > 1)
  {
    cout << __PRETTY_FUNCTION__ << " : packet" << m_packet_id << endl
         << "   data_length = " << data_length << endl
         << "   data_padding = " << data_padding << endl
         << "   dma_words_buffer = " << dma_words_buffer << endl
         << "   l2 = " << l2 << endl
         << "   dma_words = " << dma_words << endl;
  }

  // demultiplexer
  for (size_t index = 0; index < dma_words; ++index)
  {
    const dma_word& dma_word_data = buffer[index];

    if (m_verbosity > 2)
    {
      cout << __PRETTY_FUNCTION__ << " : processing DMA word "
           << index << "/" << dma_words << " with header 0x"
           << hex << dma_word_data.dma_header << dec << endl;
    }

    if ((dma_word_data.dma_header & 0xFF00) == FEE_MAGIC_KEY)
    {
      unsigned int fee_id = dma_word_data.dma_header & 0xff;

      if (fee_id < MAX_FEECOUNT)
      {
        for (unsigned int i = 0; i < DAM_DMA_WORD_LENGTH - 1; i++)
        {
          m_feeData[fee_id].push_back(dma_word_data.data[i]);
        }
        h_norm->Fill("DMA_WORD_FEE", 1);

        // immediate fee buffer processing to reduce memory consuption
        process_fee_data(fee_id);
      }
      else
      {
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE ID " << fee_id << " at position " << index << endl;
        index += DAM_DMA_WORD_LENGTH - 1;
        h_norm->Fill("DMA_WORD_FEE_INVALID", 1);
      }
    }

    else if ((dma_word_data.dma_header & 0xFF00) == GTM_MAGIC_KEY)
    {
      decode_gtm_data(dma_word_data);
      h_norm->Fill("DMA_WORD_GTM", 1);
    }
    else
    {
      cout << __PRETTY_FUNCTION__ << " : Error : Unknown data type at position " << index << ": "
           << hex << buffer[index].dma_header << dec << endl;
      // not FEE data, e.g. GTM data or other stream, to be decoded
      h_norm->Fill("DMA_WORD_INVALID", 1);
    }
  }

  // sanity check for the timeframe size
  for (auto& timeframe : m_timeFrameMap)
  {
    if (timeframe.second.size() > kMaxRawHitLimit)
    {
      cout << __PRETTY_FUNCTION__ << " : Warning : impossible amount of hits in the same timeframe at BCO "
           << timeframe.first << " : " << timeframe.second.size() << ", limit is " << kMaxRawHitLimit
           << ". Dropping this time frame!"
           << endl;
      h_norm->Fill("TimeFrameSizeLimitError", 1);

      while (!timeframe.second.empty())
      {
        delete timeframe.second.back();
        timeframe.second.pop_back();
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTimeFrameBuilder::process_fee_data(unsigned int fee)
{
  assert(m_hFEEDataStream);

  if (m_verbosity > 1)
  {
    cout << __PRETTY_FUNCTION__ << " : processing FEE " << fee << " with " << m_feeData[fee].size() << " words" << endl;
  }

  assert(fee < m_feeData.size());
  auto& data_buffer = m_feeData[fee];

  while (HEADER_LENGTH <= data_buffer.size())
  {
    // packet loop
    if (data_buffer[1] != FEE_PACKET_MAGIC_KEY_1)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE magic key at position 1 0x" << hex << data_buffer[1] << dec << endl;
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
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE magic key at position 2 0x" << hex << data_buffer[2] << dec << endl;
      }
      m_hFEEDataStream->Fill(fee, "WordSkipped", 1);
      data_buffer.pop_front();
      continue;
    }
    assert(data_buffer[2] == FEE_PACKET_MAGIC_KEY_2);

    // valid packet
    const uint16_t& pkt_length = data_buffer[0];  // this is indeed the number of 10-bit words + 5 in this packet
    if (pkt_length > MAX_PACKET_LENGTH)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE pkt_length " << pkt_length << endl;
      }
      m_hFEEDataStream->Fill(fee, "InvalidLength", 1);
      data_buffer.pop_front();
      continue;
    }

    if (pkt_length + 1U > data_buffer.size())
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << " : packet over buffer boundary for now, skip decoding and wait for more data: "
                                       " pkt_length = "
             << pkt_length
             << " data_buffer.size() = " << data_buffer.size()
             << endl;
      }
      break;
    }

    fee_payload payload;
    // continue the decoding
    payload.adc_length = data_buffer[0] - HEADER_LENGTH;  // this is indeed the number of 10-bit words in this packet
    payload.data_parity = data_buffer[4] >> 9;
    payload.sampa_address = (data_buffer[4] >> 5) & 0xf;
    payload.sampa_channel = data_buffer[4] & 0x1f;
    payload.channel = data_buffer[4] & 0x1ff;
    payload.type = (data_buffer[3] >> 7) & 0x7;
    payload.user_word = data_buffer[3] & 0x7f;
    payload.bx_timestamp = ((data_buffer[6] & 0x3ff) << 10) | (data_buffer[5] & 0x3ff);

    payload.data_crc = data_buffer[pkt_length];
    payload.calc_crc = crc16(fee, 0, pkt_length);
    if (payload.data_crc != payload.calc_crc)
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << " : CRC error in FEE "
             << fee << " at position " << pkt_length - 1
             << ": data_crc = " << payload.data_crc
             << " calc_crc = " << payload.calc_crc << endl;
      }
      m_hFEEDataStream->Fill(fee, "CRCError", 1);
      // continue;
    }
    assert(fee < m_bcoMatchingInformation_vec.size());
    auto& m_bcoMatchingInformation = m_bcoMatchingInformation_vec[fee];

    // gtm_bco matching
    if (payload.type == m_bcoMatchingInformation.HEARTBEAT_T)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__
             << " : received heartbeat packet from FEE " << fee << endl;
      }

      // if bco matching information is still not verified, drop the packet
      if (not m_bcoMatchingInformation.is_verified())
      {
        m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncUnavailable", 1);

        std::cout << "TpcTimeFrameBuilder::process_fee_data - bco_matching not verified for heart beat, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
      else  //       if (not m_bcoMatchingInformation.is_verified())
      {
        const auto result = m_bcoMatchingInformation.find_reference_heartbeat(payload);
        m_hFEEDataStream->Fill(fee, "PacketHeartBeat", 1);
        assert(m_hFEEClockAdjustment);
        m_hFEEClockAdjustment->Fill(fee, m_bcoMatchingInformation.get_multiplier_adjustment(), 1);

        if (result)
        {
          // assign gtm bco
          payload.gtm_bco = result.value();
          m_hFEEDataStream->Fill(fee, "PacketHeartBeatClockSyncOK", 1);
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
    else  //     if (payload.type == m_bcoMatchingInformation.HEARTBEAT_T)
    {
      m_hFEEChannelPacketCount->Fill(fee * MAX_CHANNELS + payload.channel, 1);

      // if bco matching information is still not verified, drop the packet
      if (not m_bcoMatchingInformation.is_verified())
      {
        m_hFEEDataStream->Fill(fee, "PacketClockSyncUnavailable", 1);

        std::cout << "TpcTimeFrameBuilder::process_fee_data - bco_matching not verified, dropping packet" << std::endl;
        m_bcoMatchingInformation.print_gtm_bco_information();
      }
      else
      {
        const auto result = m_bcoMatchingInformation.find_gtm_bco(payload.bx_timestamp);

        if (result)
        {
          // assign gtm bco
          payload.gtm_bco = result.value();
          m_hFEEDataStream->Fill(fee, "PacketClockSyncOK", 1);
        }
        else
        {
          m_hFEEDataStream->Fill(fee, "PacketClockSyncError", 1);

          // skip the waverform
        }
      }
    }

    if (m_verbosity > 2)
    {
      cout << __PRETTY_FUNCTION__ << " : received data packet "
           << " pkt_length = " << pkt_length << endl
           << " adc_length = " << payload.adc_length << endl
           << " sampa_address = " << payload.sampa_address << endl
           << " sampa_channel = " << payload.sampa_channel << endl
           << " channel = " << payload.channel << endl
           << " bx_timestamp = 0x" << hex << payload.bx_timestamp << dec << endl
           << " bco = 0x" << hex << payload.gtm_bco << dec << endl
           << " data_crc = 0x" << hex << payload.data_crc << dec << endl
           << " calc_crc = 0x" << hex << payload.calc_crc << dec << endl;
    }

    // // valid packet in the buffer, create a new hit
    // TpcRawHit *hit = new TpcRawHitv2();
    // m_timeFrameMap[gtm_bco].push_back(hit);

    // hit->set_bco(bx_timestamp);
    // hit->set_gtm_bco(gtm_bco);
    // hit->set_packetid(m_packet_id);
    // hit->set_fee(fee);
    // hit->set_channel(channel);
    // hit->set_sampaaddress(sampa_address);
    // hit->set_sampachannel(sampa_channel);
    // m_hFEEDataStream->Fill(fee, "RawHit", 1);

    // Format is (N sample) (start time), (1st sample)... (Nth sample)
    size_t pos = HEADER_LENGTH;
    while (pos < pkt_length)
    {
      const uint16_t& nsamp = data_buffer[pos++];
      const uint16_t& start_t = data_buffer[pos++];

      if (pos + nsamp >= pkt_length)
      {
        if (m_verbosity > 1)
        {
          cout << __PRETTY_FUNCTION__ << ": WARNING : nsamp: " << nsamp
               << ", pos: " << pos
               << " > pkt_length: " << pkt_length << ", format error" << endl;
        }
        m_hFEEDataStream->Fill(fee, "HitFormatError", 1);

        break;
      }

      const unsigned int fee_sampa_address = fee * MAX_SAMPA + payload.sampa_address;
      std::vector<uint16_t> adc(nsamp);
      for (int j = 0; j < nsamp; j++)
      {
        adc[j] = data_buffer[pos++];

        m_hFEESAMPAADC->Fill(start_t + j, fee_sampa_address, adc[j]);
      }
      payload.waveforms.push_back(std::make_pair(start_t, std::move(adc)));

      //   // an exception to deal with the last sample that is missing in the current hit format
      //   if (pos + 1 == pkt_length) break;
    }

    data_buffer.erase(data_buffer.begin(), data_buffer.begin() + pkt_length + 1);
    m_hFEEDataStream->Fill(fee, "WordValid", pkt_length + 1);

  }  //     while (HEADER_LENGTH < data_buffer.size())

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTimeFrameBuilder::decode_gtm_data(const TpcTimeFrameBuilder::dma_word& gtm_word)
{
  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << " : processing GTM data " << endl;
  }

  const unsigned char* gtm = reinterpret_cast<const unsigned char*>(&gtm_word);

  gtm_payload payload;

  payload.pkt_type = gtm[0] | ((unsigned short) gtm[1] << 8);
  //    if (payload.pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY && payload.pkt_type != GTM_ENDAT_MAGIC_KEY)
  if (payload.pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY && payload.pkt_type != GTM_ENDAT_MAGIC_KEY && payload.pkt_type != GTM_MODEBIT_MAGIC_KEY)
  {
    return -1;
  }

  payload.is_lvl1 = payload.pkt_type == GTM_LVL1_ACCEPT_MAGIC_KEY;
  payload.is_endat = payload.pkt_type == GTM_ENDAT_MAGIC_KEY;
  payload.is_modebit = payload.pkt_type == GTM_MODEBIT_MAGIC_KEY;

  payload.bco = ((unsigned long long) gtm[2] << 0) | ((unsigned long long) gtm[3] << 8) | ((unsigned long long) gtm[4] << 16) | ((unsigned long long) gtm[5] << 24) | ((unsigned long long) gtm[6] << 32) | (((unsigned long long) gtm[7]) << 40);
  payload.lvl1_count = ((unsigned int) gtm[8] << 0) | ((unsigned int) gtm[9] << 8) | ((unsigned int) gtm[10] << 16) | ((unsigned int) gtm[11] << 24);
  payload.endat_count = ((unsigned int) gtm[12] << 0) | ((unsigned int) gtm[13] << 8) | ((unsigned int) gtm[14] << 16) | ((unsigned int) gtm[15] << 24);
  payload.last_bco = ((unsigned long long) gtm[16] << 0) | ((unsigned long long) gtm[17] << 8) | ((unsigned long long) gtm[18] << 16) | ((unsigned long long) gtm[19] << 24) | ((unsigned long long) gtm[20] << 32) | (((unsigned long long) gtm[21]) << 40);
  payload.modebits = gtm[22];
  payload.userbits = gtm[23];

  // constexpr uint64_t bco_limit = 1ULL << GTMBCObits;
  // assert(payload.bco > bco_limit);
  // if (payload.bco < m_GTMBCOLastReading)
  // {
  //   cout << __PRETTY_FUNCTION__ << " : Info : GTM BCO rollover detected, last reading " << m_GTMBCOLastReading
  //        << " current reading " << payload.bco << " on packet " << m_packet_id << endl;
  //   ++m_GTMBCORollOverCounter;
  // }
  // m_GTMBCOLastReading = payload.bco;
  // uint64_t rollover_corrected_bco = (m_GTMBCORollOverCounter << GTMBCObits) + payload.bco;
  // m_gtmData[rollover_corrected_bco] = payload;

  if (m_verbosity > 2)
  {
    cout << " GTM data : "
         << " pkt_type = " << payload.pkt_type << endl
         << " is_lvl1 = " << payload.is_lvl1 << endl
         << " is_endat = " << payload.is_endat << endl
         << " is_modebit = " << payload.is_modebit << endl
         << " bco = 0x" << hex << payload.bco << dec << endl
         << " lvl1_count = " << payload.lvl1_count << endl
         << " endat_count = " << payload.endat_count << endl
         << " last_bco = 0x" << hex << payload.last_bco << dec << endl
         << " modebits =  0x" << hex << (int) payload.modebits << dec << endl
         << " userbits =  0x" << hex << (int) payload.userbits << dec << endl;
  }

  int fee = -1;
  for (auto& bcoMatchingInformation : m_bcoMatchingInformation_vec)
  {
    ++fee;

    if (m_verbosity > 2)
    {
      cout << __PRETTY_FUNCTION__ << " : processing GTM data for FEE " << fee << endl;
    }

    bcoMatchingInformation.save_gtm_bco_information(payload);

    if (m_verbosity > 2)
    {
      bcoMatchingInformation.print_gtm_bco_information();
    }
  }
  return 0;
}

uint16_t TpcTimeFrameBuilder::reverseBits(const uint16_t x) const
{
  uint16_t n = x;
  n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
  n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
  n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
  n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
  // n = (n >> 16) & 0x0000ffff | (n << 16) & 0xffff0000;
  return n;
}

uint16_t TpcTimeFrameBuilder::crc16(
    const unsigned int fee, const unsigned int index, const int l) const
{
  uint16_t crc = 0xffff;

  for (int i = 0; i < l; i++)
  {
    uint16_t x = m_feeData[fee][index + i];
    crc ^= reverseBits(x);
    for (uint16_t k = 0; k < 16; k++)
    {
      crc = crc & 1 ? (crc >> 1) ^ 0xa001 : crc >> 1;
    }
  }
  crc = reverseBits(crc);
  return crc;
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
      o << " }";
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
      o << " }";
    }
    return o;
  }

}  // namespace

//___________________________________________________
std::optional<uint32_t> TpcTimeFrameBuilder::BcoMatchingInformation::get_predicted_fee_bco(uint64_t gtm_bco) const
{
  // check proper initialization
  if (!is_verified())
  {
    return std::nullopt;
  }

  // get gtm bco difference with proper rollover accounting
  const uint64_t gtm_bco_difference = (gtm_bco >= m_gtm_bco_first) ? (gtm_bco - m_gtm_bco_first) : (gtm_bco + (1ULL << 40U) - m_gtm_bco_first);

  // convert to fee bco, and truncate to 20 bits
  const uint64_t fee_bco_predicted = m_fee_bco_first + get_adjusted_multiplier() * gtm_bco_difference;
  return uint32_t(fee_bco_predicted & 0xFFFFFU);
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information() const
{
  if (!m_gtm_bco_trig_list.empty())
  {
    std::cout
        << "TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information -"
        << " m_gtm_bco_trig_list: " << std::hex << m_gtm_bco_trig_list << std::dec
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
          << "TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information -"
          << " m_gtm_bco_trig_list fee predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
    }
  }

  if (!m_gtm_bco_heartbeat_list.empty())
  {
    std::cout
        << "TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information -"
        << " m_gtm_bco_heartbeat_list: " << std::hex << m_gtm_bco_heartbeat_list << std::dec
        << std::endl;

    // also print predicted fee bco
    if (is_verified())
    {
      std::list<uint32_t> fee_bco_predicted_list;
      std::transform(
          m_gtm_bco_heartbeat_list.begin(),
          m_gtm_bco_heartbeat_list.end(),
          std::back_inserter(fee_bco_predicted_list),
          [this](const uint64_t& gtm_bco) { return get_predicted_fee_bco(gtm_bco).value(); });

      std::cout
          << "TpcTimeFrameBuilder::BcoMatchingInformation::print_gtm_bco_information -"
          << " m_gtm_bco_heartbeat_list fee predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
    }
  }
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::save_gtm_bco_information(const TpcTimeFrameBuilder::gtm_payload& gtm_tagger)
{
  // append gtm_bco from taggers in this event to packet-specific list of available lv1_bco

  // save level1 trigger bco
  const bool is_lvl1 = gtm_tagger.is_lvl1;
  const bool is_endat = gtm_tagger.is_endat;
  const bool is_modebit = gtm_tagger.is_modebit;
  if (is_lvl1)
  {
    const uint64_t& gtm_bco = gtm_tagger.bco;
    m_gtm_bco_trig_list.push_back(gtm_bco);
  }

  // also save ENDDAT bco
  else if (is_endat)
  {
    const uint64_t& gtm_bco = gtm_tagger.bco;

    // add to list if difference to last entry is big enough
    if (m_gtm_bco_trig_list.empty() || (gtm_bco - m_gtm_bco_trig_list.back()) > m_max_lv1_endat_bco_diff)
    {
      m_gtm_bco_trig_list.push_back(gtm_bco);
    }
  }

  // also save hearbeat bco
  else if (is_modebit)
  {
    // get modebits
    const uint64_t& modebits = gtm_tagger.modebits;
    if (modebits & (1U << ELINK_HEARTBEAT_T))
    {
      const uint64_t& gtm_bco = gtm_tagger.bco;
      m_gtm_bco_heartbeat_list.push_back(gtm_bco);
    }

    if (modebits & (1U << BX_COUNTER_SYNC_T))
    {
      // get BCO and assign
      const uint64_t& gtm_bco = gtm_tagger.bco;
      m_gtm_bco_first = gtm_bco;
      m_fee_bco_first = 0;
      m_verified_from_modebits = true;

      if (m_verbosity)
      {
        std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::find_reference_from_modebits"
                  << " found reference from modebits BX_COUNTER_SYNC_T "
                  << "at gtm_bco = 0x" << hex << gtm_bco << dec
                  << std::endl;
      }
    }
  }
}

//___________________________________________________
std::optional<uint64_t> TpcTimeFrameBuilder::BcoMatchingInformation::find_reference_heartbeat(const TpcTimeFrameBuilder::fee_payload& HeartBeatPacket)
{
  // make sure the bco matching is properly initialized and historical valid
  if (!is_verified())
  {
    return std::nullopt;
  }

  const auto& fee_bco = HeartBeatPacket.bx_timestamp;

  // find element for which predicted fee_bco matches fee_bco, within limit
  const auto iter = std::find_if(
      m_gtm_bco_heartbeat_list.begin(),
      m_gtm_bco_heartbeat_list.end(),
      [this, fee_bco](const uint64_t& gtm_bco) {
        return get_bco_diff(get_predicted_fee_bco(gtm_bco).value(), fee_bco) < m_max_gtm_bco_diff;
      });

  // check
  if (iter != m_gtm_bco_heartbeat_list.end())
  {
    const auto gtm_bco = *iter;
    if (verbosity())
    {
      const auto fee_bco_predicted = get_predicted_fee_bco(gtm_bco).value();
      const auto fee_bco_diff = get_bco_diff(fee_bco_predicted, fee_bco);

      std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::find_reference_heartbeat -"
                << std::hex
                << " fee_bco: 0x" << fee_bco << std::endl
                << " predicted: 0x" << fee_bco_predicted << std::endl
                << " gtm_bco: 0x" << gtm_bco
                << std::dec << std::endl
                << " difference: " << fee_bco_diff << std::endl;
    }

    // update clock adjustment
    update_multiplier_adjustment(gtm_bco, fee_bco);

    return gtm_bco;
  }

  if (verbosity() > 1)
  {
    std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::find_reference_heartbeat - WARNING: failed match for fee_bco = 0x" << hex << fee_bco << dec << std::endl;
  }
  return std::nullopt;
}

//___________________________________________________
std::optional<uint64_t> TpcTimeFrameBuilder::BcoMatchingInformation::find_gtm_bco(uint32_t fee_bco)
{
  // make sure the bco matching is properly initialized
  if (!is_verified())
  {
    return std::nullopt;
  }
  // find matching gtm bco in map
  const auto bco_matching_iter = std::find_if(
      m_bco_matching_list.begin(),
      m_bco_matching_list.end(),
      [fee_bco](const m_bco_matching_pair_t& pair) { return get_bco_diff(pair.first, fee_bco) < m_max_fee_bco_diff; });

  if (bco_matching_iter != m_bco_matching_list.end())
  {
    return bco_matching_iter->second;
  }
  else
  {
    // find element for which predicted fee_bco matches fee_bco, within limit
    const auto iter = std::find_if(
        m_gtm_bco_trig_list.begin(),
        m_gtm_bco_trig_list.end(),
        [this, fee_bco](const uint64_t& gtm_bco) { return get_bco_diff(get_predicted_fee_bco(gtm_bco).value(), fee_bco) < m_max_gtm_bco_diff; });

    // check
    if (iter != m_gtm_bco_trig_list.end())
    {
      const auto gtm_bco = *iter;
      if (verbosity())
      {
        const auto fee_bco_predicted = get_predicted_fee_bco(gtm_bco).value();
        const auto fee_bco_diff = get_bco_diff(fee_bco_predicted, fee_bco);

        std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::find_gtm_bco -"
                  << std::hex
                  << " fee_bco: 0x" << fee_bco
                  << " predicted: 0x" << fee_bco_predicted
                  << " gtm_bco: 0x" << gtm_bco
                  << std::dec
                  << " difference: " << fee_bco_diff
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
      if (m_orphans.insert(fee_bco).second)
      {
        if (verbosity())
        {
          // find element for which predicted fee_bco is the closest to request
          const auto iter2 = std::min_element(
              m_gtm_bco_trig_list.begin(),
              m_gtm_bco_trig_list.end(),
              [this, fee_bco](const uint64_t& first, const uint64_t& second) { return get_bco_diff(get_predicted_fee_bco(first).value(), fee_bco) < get_bco_diff(get_predicted_fee_bco(second).value(), fee_bco); });

          const int fee_bco_diff = (iter2 != m_gtm_bco_trig_list.end()) ? get_bco_diff(get_predicted_fee_bco(*iter2).value(), fee_bco) : -1;

          std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::find_gtm_bco -"
                    << std::hex
                    << " fee_bco: 0x" << fee_bco
                    << std::dec
                    << " gtm_bco: none"
                    << " difference: " << fee_bco_diff
                    << std::endl;
        }
      }
      return std::nullopt;
    }
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
  while (m_gtm_bco_heartbeat_list.size() > m_max_matching_data_size)
  {
    m_gtm_bco_heartbeat_list.pop_front();
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
  m_gtm_bco_heartbeat_list.erase(std::remove_if(m_gtm_bco_heartbeat_list.begin(), m_gtm_bco_heartbeat_list.end(),
                                                [ref_bco](const uint64_t& bco) { return bco <= ref_bco; }),
                                 m_gtm_bco_heartbeat_list.end());

  // erase all elements from bco_list that are less than or equal to ref_bco
  m_bco_matching_list.erase(std::remove_if(m_bco_matching_list.begin(), m_bco_matching_list.end(), 
  [ref_bco](const m_bco_matching_pair_t& pair) { 
    return pair.second <= ref_bco;
     }), m_bco_matching_list.end()
     );

  // clear orphans
  m_orphans.clear();
}

//___________________________________________________
double TpcTimeFrameBuilder::BcoMatchingInformation::get_adjusted_multiplier() const
{
  return m_multiplier + m_multiplier_adjustment;
}

//___________________________________________________
void TpcTimeFrameBuilder::BcoMatchingInformation::update_multiplier_adjustment(uint64_t gtm_bco, uint32_t fee_bco)
{
  // check that references are valid
  if (!is_verified())
  {
    return;
  }

  // skip if trivial
  if (gtm_bco == m_gtm_bco_first)
  {
    return;
  }

  const uint32_t fee_bco_predicted = get_predicted_fee_bco(gtm_bco).value();
  const double delta_fee_bco = double(fee_bco) - double(fee_bco_predicted);
  const double gtm_bco_difference = (gtm_bco >= m_gtm_bco_first) ? (gtm_bco - m_gtm_bco_first) : (gtm_bco + (1ULL << 40U) - m_gtm_bco_first);

  m_multiplier_adjustment_numerator += gtm_bco_difference * delta_fee_bco;
  m_multiplier_adjustment_denominator += gtm_bco_difference * gtm_bco_difference;
  ++m_multiplier_adjustment_count;

  if (verbosity())
  {
    const auto default_precision{std::cout.precision()};
    std::cout << "TpcTimeFrameBuilder::BcoMatchingInformation::update_multiplier_adjustment -"
              << " m_multiplier_adjustment_count: " << m_multiplier_adjustment_count << endl
              << std::setprecision(10)
              << " m_multiplier: " << get_adjusted_multiplier()<< endl
              << " adjustment: " << m_multiplier_adjustment_numerator / m_multiplier_adjustment_denominator<< endl
              << " m_multiplier_adjusted: " << get_adjusted_multiplier() + m_multiplier_adjustment_numerator / m_multiplier_adjustment_denominator
              << std::setprecision(default_precision)
              << std::endl;
  }

  // update multiplier
  if (m_multiplier_adjustment_count > m_max_multiplier_adjustment_count)
  {
    m_multiplier_adjustment += m_multiplier_adjustment_numerator / m_multiplier_adjustment_denominator;
    m_multiplier_adjustment_numerator = 0;
    m_multiplier_adjustment_denominator = 0;
    m_multiplier_adjustment_count = 0;
  }
}

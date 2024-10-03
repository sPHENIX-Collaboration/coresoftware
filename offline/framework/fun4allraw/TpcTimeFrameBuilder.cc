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
{
  m_feeData.resize(MAX_FEECOUNT);

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h = new TH1D(TString(m_HistoPrefix.c_str()) + "_Normalization",  //
                    TString(m_HistoPrefix.c_str()) + " Normalization;Items;Count", 10, .5, 10.5);
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

  h->GetYaxis()->SetBinLabel(i++, "TimeFrameSizeLimitError");
  assert(i <= 10);
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
                              MAX_FEECOUNT, -.5, MAX_FEECOUNT - .5, 2, .5, 10.5);
  i = 1;
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordValid");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "WordSkipped");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "RawHit");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitFormatError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "MissingLastHit");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitCRCError");
  m_hFEEDataStream->GetYaxis()->SetBinLabel(i++, "HitUnusedBeforeCleanup");
  assert(i <= 10);
  hm->registerHisto(m_hFEEDataStream);
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

bool TpcTimeFrameBuilder::isMoreDataRequired() const
{
  if (m_gtmData.size() == 0) return true;

  if (m_gtmData.rbegin()->first - m_gtmData.begin()->first > kFEEDataTransmissionWindow)
    return false;

  return true;
}

void TpcTimeFrameBuilder::CleanupUsedPackets(const uint64_t bclk)
{
  if (m_verbosity > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << " packet " << m_packet_id << ": cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }

  uint64_t bclk_rollover_corrected = bclk;
  if (m_gtmData.begin() != m_gtmData.end())
  {
    if ((m_gtmData.begin()->first & GTMBCOmask_ValidBits) > bclk + (1ULL << (GTMBCObits - 1)))
    {
      bclk_rollover_corrected = (((m_gtmData.begin()->first >> GTMBCObits) + 1) << GTMBCObits) | (bclk & GTMBCOmask_ValidBits);
    }
    else if ((m_gtmData.begin()->first & GTMBCOmask_ValidBits) + (1ULL << (GTMBCObits - 1)) < bclk)
    {
      bclk_rollover_corrected = (((m_gtmData.begin()->first >> GTMBCObits) - 1) << GTMBCObits) | (bclk & GTMBCOmask_ValidBits);
    }
    else
    {
      bclk_rollover_corrected = (m_gtmData.begin()->first & GTMBCOmask_RollOverCounts) | (bclk & GTMBCOmask_ValidBits);
    }
  }

  assert(m_hFEEDataStream);

  for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)
  {
    if (it->first <= bclk_rollover_corrected)
    {
      while (!it->second.empty())
      {
        m_hFEEDataStream->Fill(it->second.back()->get_fee(), "HitUnusedBeforeCleanup", 1);
        delete it->second.back();
        it->second.pop_back();
      }
      m_timeFrameMap.erase(it++);
    }
    else
    {
      break;
    }
  }  //   for (auto it = m_timeFrameMap.begin(); it != m_timeFrameMap.end();)

  for (auto it = m_gtmData.begin(); it != m_gtmData.end();)
  {
    if (it->first <= bclk)
    {
      m_gtmData.erase(it++);
    }
    else
    {
      break;
    }
  }  //   for (auto it = m_gtmData.begin(); it != m_gtmData.end();)
}

int TpcTimeFrameBuilder::ProcessPacket(Packet *packet)
{
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

  if (m_verbosity)
  {
    cout << __PRETTY_FUNCTION__ << " : received packet ";
    packet->identify();
  }

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
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
  packet->fillIntArray(reinterpret_cast<int *>(buffer.data()), data_length + DAM_DMA_WORD_LENGTH / 2, &l2, "DATA");
  assert(l2 <= data_length);
  l2 -= data_padding;
  assert(l2 >= 0);

  size_t dma_words = l2 * 2 / DAM_DMA_WORD_LENGTH;
  assert(dma_words <= buffer.size());
  assert(h_PacketLength_Residual);
  h_PacketLength_Residual->Fill((l2 * 2) % DAM_DMA_WORD_LENGTH);

  // demultiplexer
  for (size_t index = 0; index < dma_words; ++index)
  {
    const dma_word &dma_word_data = buffer[index];

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

  for (int fee_id = 0; fee_id < MAX_FEECOUNT; ++fee_id)
  {
    process_fee_data(fee_id);
  }

  // sanity check for the timeframe size
  for (auto &timeframe : m_timeFrameMap)
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

  if (m_verbosity)
  {
    cout << __PRETTY_FUNCTION__ << " : processing FEE " << fee << " with " << m_feeData[fee].size() << " words" << endl;
  }

  assert(fee < m_feeData.size());
  auto &data_buffer = m_feeData[fee];

  while (HEADER_LENGTH <= data_buffer.size())
  {
    // packet loop
    if (data_buffer[1] != FEE_PACKET_MAGIC_KEY_1)
    {
      if (m_verbosity > 1)
      {
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE magic key at position 1 " << data_buffer[1] << endl;
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
        cout << __PRETTY_FUNCTION__ << " : Error : Invalid FEE magic key at position 2 " << data_buffer[2] << endl;
      }
      m_hFEEDataStream->Fill(fee, "WordSkipped", 1);
      data_buffer.pop_front();
      continue;
    }
    assert(data_buffer[2] == FEE_PACKET_MAGIC_KEY_2);

    // valid packet
    const uint16_t &pkt_length = data_buffer[0];  // this is indeed the number of 10-bit words + 5 in this packet
    if (pkt_length > data_buffer.size())
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

    // continue the decoding
    const uint16_t adc_length = data_buffer[0] - HEADER_LENGTH;  // this is indeed the number of 10-bit words in this packet
    // const uint16_t data_parity = data_buffer[4] >> 9;
    const uint16_t sampa_address = (data_buffer[4] >> 5) & 0xf;
    const uint16_t sampa_channel = data_buffer[4] & 0x1f;
    const uint16_t channel = data_buffer[4] & 0x1ff;
    // const uint16_t type = (data_buffer[3] >> 7) & 0x7;
    // const uint16_t user_word = data_buffer[3] & 0x7f;
    const uint16_t bx_timestamp = ((data_buffer[6] & 0x3ff) << 10) | (data_buffer[5] & 0x3ff);

    if (m_verbosity > 1)
    {
      cout << __PRETTY_FUNCTION__ << " : received data packet "
           << " pkt_length = " << pkt_length
           << " adc_length = " << adc_length
           << " sampa_address = " << sampa_address
           << " sampa_channel = " << sampa_channel
           << " channel = " << channel
           << " bx_timestamp = " << bx_timestamp

           << endl;
    }

    // gtm_bco matching
    uint64_t gtm_bco = matchFEE2GTMBCO(bx_timestamp);

    // valid packet in the buffer, create a new hit
    TpcRawHit *hit = new TpcRawHitv2();
    m_timeFrameMap[gtm_bco].push_back(hit);

    hit->set_bco(bx_timestamp);
    hit->set_gtm_bco(gtm_bco);
    hit->set_packetid(m_packet_id);
    hit->set_fee(fee);
    hit->set_channel(channel);
    hit->set_sampaaddress(sampa_address);
    hit->set_sampachannel(sampa_channel);
    m_hFEEDataStream->Fill(fee, "RawHit", 1);

    const uint16_t data_crc = data_buffer[pkt_length - 1];
    const uint16_t calc_crc = crc16(fee, 0, pkt_length - 1);
    if (data_crc != calc_crc)
    {
      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << " : CRC error in FEE " << fee << " at position " << pkt_length - 1 << ": data_crc = " << data_crc << " calc_crc = " << calc_crc << endl;
      }
      m_hFEEDataStream->Fill(fee, "CRCError", 1);
      continue;
    }

    size_t pos = HEADER_LENGTH;
    // Format is (N sample) (start time), (1st sample)... (Nth sample)
    while (pos < pkt_length)
    {
      const uint16_t &nsamp = data_buffer[pos++];
      const uint16_t &start_t = data_buffer[pos++];

      if (pos + nsamp >= pkt_length)
      {
        if (m_verbosity > 1)
        {
          cout << __PRETTY_FUNCTION__ << ": WARNING : nsamp: " << nsamp
               << ", pos: " << pos
               << " > pkt_length: " << pkt_length << ", format error" << endl;
        }
        m_hFEEDataStream->Fill(fee, "HitFormatError", 1);

        delete hit;

        break;
      }

      for (int j = 0; j < nsamp; j++)
      {
        hit->set_adc(start_t + j, data_buffer[pos++]);

        // an exception to deal with the last sample that is missing in the current hit format
        if (pos + 1 == pkt_length)
        {
          m_hFEEDataStream->Fill(fee, "MissingLastADC", 1);
          break;
        }
      }

      // an exception to deal with the last sample that is missing in the current hit format
      if (pos + 1 == pkt_length) break;
    }

    data_buffer.erase(data_buffer.begin(), data_buffer.begin() + pkt_length);
    m_hFEEDataStream->Fill(fee, "WordValid", pkt_length);

  }  //     while (HEADER_LENGTH < data_buffer.size())

  return Fun4AllReturnCodes::EVENT_OK;
}

uint64_t TpcTimeFrameBuilder::matchFEE2GTMBCO(uint16_t fee_bco)
{
  if (m_verbosity > 2)
  {
    cout << __PRETTY_FUNCTION__ << " : FEE BCO " << fee_bco << endl;
  }

  uint64_t gtm_bco = 0;
  switch (m_gtmMatcherStrategy)
  {
  case kLastLv1Tagger:
  {
    // find the last GTM BCO that is less than the FEE BCO
    if (m_gtmData.rbegin() != m_gtmData.rend())
    {
      gtm_bco = m_gtmData.rbegin()->first;

      if (m_verbosity > 2)
      {
        cout << __PRETTY_FUNCTION__ << " : strategy kLastLv1Tagger, and match to gtm_bco " << gtm_bco << endl;
      }
    }
    break;
  }
  case kFEEWaveformBCOSync:
  {
    // // find the GTM BCO that is closest to the FEE BCO
    // uint64_t min_diff = 0xffffffffffffffff;
    // for (auto it = m_gtmData.begin(); it != m_gtmData.end(); ++it)
    // {
    //   uint64_t diff = it->first - fee_bco;
    //   if (diff < min_diff)
    //   {
    //     min_diff = diff;
    //     gtm_bco = it->first;
    //   }
    // }
    break;
  }
  // case kFEEHeartBeatSync:
  // {
  //   break;
  // }
  default:
  {
    cout << __PRETTY_FUNCTION__ << " : Error : Unknown GTM Matcher Strategy " << m_gtmMatcherStrategy << endl;
    assert(0);
  }
  }

  return gtm_bco;
}

int TpcTimeFrameBuilder::decode_gtm_data(const TpcTimeFrameBuilder::dma_word &gtm_word)
{
  const unsigned char *gtm = reinterpret_cast<const unsigned char *>(&gtm_word);

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

  constexpr uint64_t bco_limit = 1ULL << GTMBCObits;
  assert(payload.bco > bco_limit);
  if (payload.bco < m_GTMBCOLastReading)
  {
    cout << __PRETTY_FUNCTION__ << " : Info : GTM BCO rollover detected, last reading " << m_GTMBCOLastReading
         << " current reading " << payload.bco << " on packet " << m_packet_id << endl;
    ++m_GTMBCORollOverCounter;
  }
  m_GTMBCOLastReading = payload.bco;
  uint64_t rollover_corrected_bco = (m_GTMBCORollOverCounter << GTMBCObits) + payload.bco;
  m_gtmData[rollover_corrected_bco] = payload;

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

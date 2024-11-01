#include "SingleMicromegasPoolInput_v2.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/MicromegasRawHitContainerv1.h>
#include <ffarawobjects/MicromegasRawHitv2.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/oncsSubConstants.h>

#include <TFile.h>
#include <TH1.h>

#include <algorithm>
#include <bitset>

namespace
{
  // maximum number of packets
  static constexpr int m_npackets_active = 2;

  // minimum number of requested samples
  static constexpr int m_min_req_samples = 5;

  /* see: https://git.racf.bnl.gov/gitea/Instrumentation/sampa_data/src/branch/fmtv2/README.md */
  enum SampaDataType
  {
    HEARTBEAT_T = 0b000,
    TRUNCATED_DATA_T = 0b001,
    TRUNCATED_TRIG_EARLY_DATA_T = 0b011,
    NORMAL_DATA_T = 0b100,
    LARGE_DATA_T = 0b101,
    TRIG_EARLY_DATA_T = 0b110,
    TRIG_EARLY_LARGE_DATA_T = 0b111,
  };

  //_____________________________________________________________
  [[maybe_unused]] uint16_t reverseBits(const uint16_t& x)
  {
    uint16_t n = x;
    n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
    n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
    n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
    n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
    // n = (n >> 16) & 0x0000ffff | (n << 16) & 0xffff0000;
    return n;
  }

  //_____________________________________________________________
  [[maybe_unused]] uint16_t crc16( const std::deque<uint16_t>& data, const unsigned int index, const int l)
  {
    uint16_t crc = 0xffff;

    for (int i = 0; i < l; i++)
    {
      uint16_t x = data[index + i];
      crc ^= reverseBits(x);
      for (uint16_t k = 0; k < 16; k++)
      {
        crc = crc & 1 ? (crc >> 1) ^ 0xa001 : crc >> 1;
      }
    }
    crc = reverseBits(crc);
    return crc;
  }

}  // namespace

//______________________________________________________________
SingleMicromegasPoolInput_v2::SingleMicromegasPoolInput_v2(const std::string& name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::MICROMEGAS);
  m_rawHitContainerName = "MICROMEGASRAWHIT";
}

//______________________________________________________________
SingleMicromegasPoolInput_v2::~SingleMicromegasPoolInput_v2()
{
  std::cout << "SingleMicromegasPoolInput_v2::~SingleMicromegasPoolInput_v2 - runnumber: " << RunNumber() << std::endl;

  // timer statistics
  m_timer.print_stat();

  // dropped waveforms
  for( const auto& [packet,counts]:m_waveform_count_total )
  {
    const auto dropped_bco =  m_waveform_count_dropped_bco[packet];
    const auto dropped_pool =  m_waveform_count_dropped_pool[packet];
    std::cout << "SingleMicromegasPoolInput_v2::~SingleMicromegasPoolInput_v2 -"
      << " packet: " << packet
      << " wf_total: " << counts
      << " wf_dropped_bco: " << dropped_bco
      << " wf_dropped_pool: " << dropped_pool
      << " ratio_bco: " << double(dropped_bco)/counts
      << " ratio_pool: " << double(dropped_pool)/counts
      << std::endl;
  }
  std::cout << std::endl;

  // drop per fee statistics
  for( const auto& [fee,counts]:m_fee_waveform_count_total )
  {
    const auto dropped_bco =  m_fee_waveform_count_dropped_bco[fee];
    std::cout << "SingleMicromegasPoolInput_v2::~SingleMicromegasPoolInput_v2 -"
      << " fee: " << fee
      << " wf_total: " << counts
      << " wf_dropped_bco: " << dropped_bco
      << " ratio_bco: " << double(dropped_bco)/counts
      << std::endl;
  }
  std::cout << std::endl;

  // also printout adjusted multipliers
  for( const auto& [packet_id, bco_matching]:m_bco_matching_information_map )
  {
    const auto old_precision = std::cout.precision();
    std::cout << "SingleMicromegasPoolInput_v2::~SingleMicromegasPoolInput_v2 -"
      << " packet: " << packet_id
      << " multiplier: " <<  std::setprecision(9) << bco_matching.get_adjusted_multiplier() << std::setprecision(old_precision)
      << std::endl;
  }

}

//______________________________________________________________
void SingleMicromegasPoolInput_v2::FillPool(const unsigned int /*nbclks*/)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }

  while (!GetEventiterator())  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }

  while (GetSomeMoreEvents())
  {
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }

      // get next event
      evt.reset(GetEventiterator()->getNextEvent());
    }

    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }

    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }

    RunNumber(evt->getRunNumber());
    if (Verbosity() > 1)
    {
      evt->identify();
    }

    const int npackets = evt->getPacketList(&plist[0], 10);

    if (npackets == 10)
    {
      std::cout << "SingleMicromegasPoolInput_v2::FillPool - too many packets" << std::endl;
      exit(1);
    }

    m_timer.restart();
    for (int i = 0; i < npackets; i++)
    {
      // keep pointer to local packet
      std::unique_ptr<Packet> packet(plist[i]);

      // process
      process_packet( packet.get() );
    }
    m_timer.stop();
  }
}

//______________________________________________________________
void SingleMicromegasPoolInput_v2::Print(const std::string& what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto& bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "FEM: " << feeiter << std::endl;
      }
    }
  }

  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }

  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto& bcliter : m_MicromegasRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout
            << "fee: " << feeiter->get_fee()
            << " at " << std::hex << feeiter << std::dec
            << std::endl;
      }
    }
  }

  if (what == "ALL" || what == "STACK")
  {
    for (const auto& iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

//____________________________________________________________________________
void SingleMicromegasPoolInput_v2::CleanupUsedPackets_with_qa(const uint64_t bclk, bool dropped)
{

  // delete all raw hits associated to bco smaller than reference, and remove from map
  for(auto iter = m_MicromegasRawHitMap.begin(); iter != m_MicromegasRawHitMap.end() && (iter->first <= bclk); iter = m_MicromegasRawHitMap.erase(iter))
  {
    for (const auto& rawhit : iter->second)
    {
      if( dropped )
      {
        // increment dropped waveform counter and histogram
        ++m_waveform_count_dropped_pool[rawhit->get_packetid()];
        h_waveform_count_dropped_pool->Fill( std::to_string(rawhit->get_packetid()).c_str(), 1 );
      }
      delete rawhit;
    }
  }

  // cleanup bco stacks
  /* it erases all elements for which the bco is no greater than the provided one */
  m_BclkStack.erase(m_BclkStack.begin(), m_BclkStack.upper_bound(bclk));
  m_BeamClockFEE.erase(m_BeamClockFEE.begin(), m_BeamClockFEE.upper_bound(bclk));
  m_BeamClockPacket.erase(m_BeamClockPacket.begin(), m_BeamClockPacket.upper_bound(bclk));

  // cleanup matching information
  for( auto&& bco_matching:m_bco_matching_information_map )
  { bco_matching.second.cleanup(bclk); }

}

//_______________________________________________________
void SingleMicromegasPoolInput_v2::ClearCurrentEvent()
{
  std::cout << "SingleMicromegasPoolInput_v2::ClearCurrentEvent." << std::endl;
  uint64_t currentbclk = *m_BclkStack.begin();
  CleanupUsedPackets(currentbclk);
  return;
}

//_______________________________________________________
bool SingleMicromegasPoolInput_v2::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }

  // check minimum pool size
  if (m_MicromegasRawHitMap.size() < m_BcoPoolSize)
  {
    return true;
  }

  // make sure that the latest BCO received by each FEEs is past the current BCO
  std::set<int> toerase;
  uint64_t lowest_bclk = m_MicromegasRawHitMap.begin()->first + m_BcoRange;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= lowest_bclk)
    {
      uint64_t highest_bclk = m_MicromegasRawHitMap.rbegin()->first;
      if ((highest_bclk - m_MicromegasRawHitMap.begin()->first) < MaxBclkDiff())
      {
        return true;
      }
      else
      {
        std::cout << PHWHERE << Name() << ": erasing FEE " << bcliter.first
                  << " with stuck bclk: " << std::hex << bcliter.second
                  << " current bco range: 0x" << m_MicromegasRawHitMap.begin()->first
                  << ", to: 0x" << highest_bclk << ", delta: " << std::dec
                  << (highest_bclk - m_MicromegasRawHitMap.begin()->first)
                  << std::dec << std::endl;
        toerase.insert(bcliter.first);
      }
    }
  }
  for(auto& fee: toerase)
  {
    m_FEEBclkMap.erase(fee);
  }

  return false;
}

//_______________________________________________________
void SingleMicromegasPoolInput_v2::CreateDSTNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  PHNodeIterator iterDst(dstNode);
  auto detNode = dynamic_cast<PHCompositeNode*>(iterDst.findFirst("PHCompositeNode", "MICROMEGAS"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("MICROMEGAS");
    dstNode->addNode(detNode);
  }

  auto container = findNode::getClass<MicromegasRawHitContainer>(detNode, m_rawHitContainerName);
  if (!container)
  {
    container = new MicromegasRawHitContainerv1();
    auto newNode = new PHIODataNode<PHObject>(container, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput_v2::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetMicromegasBcoRange(m_BcoRange);
    StreamingInputManager()->SetMicromegasNegativeBco(m_NegativeBco);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput_v2::FillBcoQA(uint64_t gtm_bco)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // set of packets for which the gtm BCO is found at least once
  unsigned int n_waveforms = 0;
  std::set<int> found_packets;
  for (uint64_t gtm_bco_loc = gtm_bco - m_NegativeBco; gtm_bco_loc < gtm_bco + m_BcoRange - m_NegativeBco; ++gtm_bco_loc)
  {
    // packet ids
    const auto packet_iter = m_BeamClockPacket.find(gtm_bco_loc);
    if( packet_iter != m_BeamClockPacket.end() )
    { found_packets.insert( packet_iter->second.begin(), packet_iter->second.end() ); }

    // waveforms
    const auto wf_iter = m_MicromegasRawHitMap.find(gtm_bco_loc);
    if (wf_iter != m_MicromegasRawHitMap.end())
    { n_waveforms += wf_iter->second.size(); }
  }

  // per packet statistics
  h_packet_stat->Fill( "Reference", 1 );

  for( const auto& packet_id:found_packets )
  { h_packet_stat->Fill( std::to_string(packet_id).c_str(), 1 ); }
  h_packet_stat->Fill( "All", found_packets.size()>= m_npackets_active );

  // how many packet_id found for this BCO
  h_packet->Fill(found_packets.size());

  // how many waveforms found for this BCO
  h_waveform->Fill(n_waveforms);
}
//_______________________________________________________
void SingleMicromegasPoolInput_v2::createQAHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // number of packets found with BCO from felix matching reference BCO
  h_packet = new TH1F( "h_MicromegasBCOQA_npacket_bco", "TPOT Packet Count per GTM BCO; Matching BCO tagger count; GL1 trigger count", 10, 0, 10 );

  // number of waveforms found with BCO from felix matching reference BCO
  h_waveform = new TH1F( "h_MicromegasBCOQA_nwaveform_bco", "TPOT Waveform Count per GTM BCO; Matching Waveform count; GL1 trigger count", 4100, 0, 4100 );

  /*
   * first bin is the number of requested GL1 BCO, for reference
   * next two bins is the number of times the GL1 BCO is found in the taggers list for a given packet_id
   * last bin is the sum
   */
  h_packet_stat = new TH1F( "h_MicromegasBCOQA_packet_stat", "Matching Tagger count per packet; packet id; GL1 trigger count", m_npackets_active+2, 0, m_npackets_active+2 );
  h_packet_stat->GetXaxis()->SetBinLabel(1, "Reference" );
  h_packet_stat->GetXaxis()->SetBinLabel(2, "5001" );
  h_packet_stat->GetXaxis()->SetBinLabel(3, "5002" );
  h_packet_stat->GetXaxis()->SetBinLabel(4, "All" );
  h_packet_stat->GetYaxis()->SetTitle( "trigger count" );

  // total number of waveform per packet
  h_waveform_count_total = new TH1F( "h_MicromegasBCOQA_waveform_count_total", "Total number of waveforms per packet", m_npackets_active, 0, m_npackets_active );

  // number of dropped waveform per packet due to bco mismatch
  h_waveform_count_dropped_bco = new TH1F( "h_MicromegasBCOQA_waveform_count_dropped_bco", "Number of dropped waveforms per packet (bco)", m_npackets_active, 0, m_npackets_active );

  // number of dropped waveform per packet due to fun4all pool mismatch
  h_waveform_count_dropped_pool = new TH1F( "h_MicromegasBCOQA_waveform_count_dropped_pool", "Number of dropped waveforms per packet (pool)", m_npackets_active, 0, m_npackets_active );

  // define axis
  for( const auto& h:std::initializer_list<TH1*>{h_waveform_count_total, h_waveform_count_dropped_bco, h_waveform_count_dropped_pool} )
  {
    h->GetXaxis()->SetBinLabel(1, "5001" );
    h->GetXaxis()->SetBinLabel(2, "5002" );
    h->GetYaxis()->SetTitle( "waveform count" );
    h->GetXaxis()->SetTitle( "packet id" );
  }

  // register all histograms to histogram manager
  for( const auto& h:std::initializer_list<TH1*>{h_packet, h_waveform, h_packet_stat, h_waveform_count_total, h_waveform_count_dropped_bco, h_waveform_count_dropped_pool} )
  {
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow);
    hm->registerHisto(h);
  }
}

//__________________________________________________________________________________
void SingleMicromegasPoolInput_v2::process_packet(Packet* packet )
{

  // check hit format
  if (packet->getHitFormat() != IDTPCFEEV4) return;

  // get packet id
  const int packet_id = packet->getIdentifier();

  // get relevant bco matching information and initialize
  auto& bco_matching_information = m_bco_matching_information_map[packet_id];
  bco_matching_information.set_verbosity(Verbosity());

  // decode
  const int data_length = packet->getDataLength();  // 32bit length
  int data_padding = packet->getPadding();  // 32bit padding

  // maximum number of dma words
  const size_t dma_words_buffer = data_length * 2 / DAM_DMA_WORD_LENGTH + 1;

  // dma words
  std::vector<dma_word> buffer(dma_words_buffer);

  int l2 = 0;
  packet->fillIntArray(reinterpret_cast<int*>(buffer.data()), data_length + DAM_DMA_WORD_LENGTH / 2, &l2, "DATA");

  // sanity checks
  assert(l2 <= data_length);

  // sanity checks
  l2 -= data_padding;
  assert(l2 >= 0);

  // actual number of dma words
  const size_t dma_words = l2 * 2 / DAM_DMA_WORD_LENGTH;
  assert(dma_words <= buffer.size());

//   // residual data (dropped)
//   const size_t dma_residual = (l2 * 2) % DAM_DMA_WORD_LENGTH;

  // demultiplexer
  for (size_t index = 0; index < dma_words; ++index)
  {
    const auto& dma_word_data = buffer[index];

    if ((dma_word_data.dma_header & 0xFF00) == GTM_MAGIC_KEY)
    {

      // decode gtm data
      decode_gtm_data(packet_id, dma_word_data);

    } else if ((dma_word_data.dma_header & 0xFF00) == FEE_MAGIC_KEY) {

      // decode fee data
      // get fee id
      const unsigned int fee_id = dma_word_data.dma_header & 0xff;

      // populate fee buffer
      if (fee_id < MAX_FEECOUNT)
      {
        for (unsigned int i = 0; i < DAM_DMA_WORD_LENGTH - 1; i++)
        { m_feeData[fee_id].push_back(dma_word_data.data[i]); }

        // immediate fee buffer processing to reduce memory consuption
        process_fee_data(packet_id, fee_id);
      }
    }
  }
}

//____________________________________________________________________
void SingleMicromegasPoolInput_v2::decode_gtm_data( int packet_id, const SingleMicromegasPoolInput_v2::dma_word& gtm_word )
{
  const unsigned char* gtm = reinterpret_cast<const unsigned char*>(&gtm_word);
  MicromegasBcoMatchingInformation_v2::gtm_payload payload;

  payload.pkt_type = gtm[0] | ((unsigned short) gtm[1] << 8);

  // check packet type
  if (payload.pkt_type != GTM_LVL1_ACCEPT_MAGIC_KEY &&
    payload.pkt_type != GTM_ENDAT_MAGIC_KEY &&
    payload.pkt_type != GTM_MODEBIT_MAGIC_KEY)
  { return; }

  payload.is_lvl1 = payload.pkt_type == GTM_LVL1_ACCEPT_MAGIC_KEY;
  payload.is_endat = payload.pkt_type == GTM_ENDAT_MAGIC_KEY;
  payload.is_modebit = payload.pkt_type == GTM_MODEBIT_MAGIC_KEY;

  payload.bco = ((unsigned long long) gtm[2] << 0) | ((unsigned long long) gtm[3] << 8) | ((unsigned long long) gtm[4] << 16) | ((unsigned long long) gtm[5] << 24) | ((unsigned long long) gtm[6] << 32) | (((unsigned long long) gtm[7]) << 40);
  payload.lvl1_count = ((unsigned int) gtm[8] << 0) | ((unsigned int) gtm[9] << 8) | ((unsigned int) gtm[10] << 16) | ((unsigned int) gtm[11] << 24);
  payload.endat_count = ((unsigned int) gtm[12] << 0) | ((unsigned int) gtm[13] << 8) | ((unsigned int) gtm[14] << 16) | ((unsigned int) gtm[15] << 24);
  payload.last_bco = ((unsigned long long) gtm[16] << 0) | ((unsigned long long) gtm[17] << 8) | ((unsigned long long) gtm[18] << 16) | ((unsigned long long) gtm[19] << 24) | ((unsigned long long) gtm[20] << 32) | (((unsigned long long) gtm[21]) << 40);
  payload.modebits = gtm[22];
  payload.userbits = gtm[23];

  // save bco information
  auto& bco_matching_information = m_bco_matching_information_map[packet_id];
  bco_matching_information.save_gtm_bco_information(payload);

  if(payload.is_lvl1||payload.is_endat)
  {
    const auto& gtm_bco = payload.bco;
    m_BeamClockPacket[gtm_bco].insert(packet_id);
    m_BclkStack.insert(gtm_bco);
  }

  // find reference from modebits, using BX_COUNTER_SYNC_T
  /*
  * This needs to be done even if the bco matching information is already verified
  * because any BX_COUNTER_SYNC_T event will break past references
  */
  bco_matching_information.find_reference_from_modebits(payload);
}

//____________________________________________________________________
void SingleMicromegasPoolInput_v2::process_fee_data( int packet_id, unsigned int fee_id )
{
  // get bco information
  auto& bco_matching_information = m_bco_matching_information_map[packet_id];

  // get fee buffer
  assert(fee_id < m_feeData.size());
  auto& data_buffer = m_feeData[fee_id];

  // process header
  // while (HEADER_LENGTH <= data_buffer.size())
  while (HEADER_LENGTH <= data_buffer.size())
  {

    // packet loop
    /* question: why check index [1] ? */
    if (data_buffer[1] != FEE_PACKET_MAGIC_KEY_1)
    {
      if (Verbosity())
      {
        std::cout << "SingleMicromegasPoolInput_v2::process_fee_data - Invalid FEE magic key at position 1 0x" << std::hex << data_buffer[1] << std::dec << std::endl;
      }
      data_buffer.pop_front();
      continue;
    }

    if (data_buffer[2] != FEE_PACKET_MAGIC_KEY_2)
    {
      if (Verbosity())
      {
        std::cout << "SingleMicromegasPoolInput_v2::process_fee_data - Invalid FEE magic key at position 2 0x" << std::hex << data_buffer[2] << std::dec << std::endl;
      }
      data_buffer.pop_front();
      continue;
    }

    // packet length
    // number of 10-bit words + 5 in this packet
    const uint16_t& pkt_length = data_buffer[0];
    if (pkt_length > MAX_PACKET_LENGTH)
    {
      if (Verbosity())
      {
        std::cout << "SingleMicromegasPoolInput_v2::process_fee_data - Invalid FEE pkt_length " << pkt_length << std::endl;
      }
      data_buffer.pop_front();
      continue;
    }

    // compare to data size
    if (pkt_length + 1U > data_buffer.size())
    {
      // not enough data. will wait for more
      break;
    }

    // increment number of waveforms
    ++m_waveform_count_total[packet_id];
    ++m_fee_waveform_count_total[fee_id];
    h_waveform_count_total->Fill( std::to_string(packet_id).c_str(), 1 );

    // create payload
    MicromegasBcoMatchingInformation_v2::fee_payload payload;

    // number of 10-bit words in this packet
    payload.adc_length = data_buffer[0] - HEADER_LENGTH;
    payload.data_parity = data_buffer[4] >> 9;
    payload.sampa_address = (data_buffer[4] >> 5) & 0xf;
    payload.sampa_channel = data_buffer[4] & 0x1f;
    payload.channel = data_buffer[4] & 0x1ff;
    payload.type = (data_buffer[3] >> 7) & 0x7;
    payload.user_word = data_buffer[3] & 0x7f;
    payload.bx_timestamp = ((data_buffer[6] & 0x3ff) << 10) | (data_buffer[5] & 0x3ff);

    // crc
    payload.data_crc = data_buffer[pkt_length];
    // payload.calc_crc = crc16(data_buffer, 0, pkt_length);

    // data
    // Format is (N sample) (start time), (1st sample)... (Nth sample)
    size_t pos = HEADER_LENGTH;
    while (pos < pkt_length)
    {
      const uint16_t& samples = data_buffer[pos++];
      const uint16_t& start_t = data_buffer[pos++];
      if (pos + samples > pkt_length)
      {
        if (Verbosity())
        {
          std::cout << "SingleMicromegasPoolInput_v2::process_fee_data -"
            << " samples: " << samples
            << " pos: " << pos
            << "pkt_length: " << pkt_length
            << " format error"
            << std::endl;
        }
        break;
      }

      std::vector<uint16_t> adc(samples);
      for (int i = 0; i < samples; ++i)
      { adc[i] = data_buffer[pos++]; }

      // add
      payload.waveforms.emplace_back(start_t,std::move(adc));
    }

    // cleanup
    data_buffer.erase(data_buffer.begin(), data_buffer.begin() + pkt_length + 1);

    // check bco matching information
    if (!bco_matching_information.is_verified())
    {
      bco_matching_information.find_reference_from_data(payload);
    }

    // if bco matching information is still not verified, drop the packet
    if (!bco_matching_information.is_verified())
    {
      ++m_waveform_count_dropped_bco[packet_id];
      h_waveform_count_dropped_bco->Fill( std::to_string(packet_id).c_str(), 1);
      continue;
    }

    // try get gtm bco matching fee
    const auto& fee_bco = payload.bx_timestamp;

    // find matching gtm bco
    uint64_t gtm_bco = 0;
    const auto result = bco_matching_information.find_gtm_bco(fee_bco);
    if (result)
    {
      // assign gtm bco
      gtm_bco = result.value();

    } else {

      // increment counter and histogram
      ++m_waveform_count_dropped_bco[packet_id];
      ++m_fee_waveform_count_dropped_bco[fee_id];
      h_waveform_count_dropped_bco->Fill( std::to_string(packet_id).c_str(), 1 );

      // skip the waverform
      continue;
    }

    // ignore heartbeat waveforms
    if( payload.type == HEARTBEAT_T )
    { continue; }

    // create new hit
    auto newhit = std::make_unique<MicromegasRawHitv2>();
    newhit->set_bco(fee_bco);
    newhit->set_gtm_bco(gtm_bco);

    // packet id, fee id, channel, etc.
    newhit->set_packetid(packet_id);
    newhit->set_fee(fee_id);
    newhit->set_channel(payload.channel);
    newhit->set_sampaaddress(payload.sampa_address);
    newhit->set_sampachannel(payload.sampa_channel);

    // assign samples
    newhit->set_sample_begin(0);
    newhit->set_sample_end(MAX_SAMPLE);

    // adc values
    for( const auto& [start_t, adc]:payload.waveforms )
    {
      for( size_t i=0; i < adc.size(); ++i )
      {
        if( start_t+i < MAX_SAMPLE )
        {
          newhit->set_adc(start_t+i,adc[i]);
        } else {
          std::cout << "SingleMicromegasPoolInput_v2::process_fee_data -"
            << " fee_id: " << fee_id
            << " channel: " << payload.channel
            << " invalid sample: " << start_t+i
            << std::endl;
          // break;
        }
      }
    }

    m_BeamClockFEE[gtm_bco].insert(fee_id);
    m_FEEBclkMap[fee_id] = gtm_bco;

    if (StreamingInputManager())
    {
      StreamingInputManager()->AddMicromegasRawHit(gtm_bco, newhit.get());
    }

    m_MicromegasRawHitMap[gtm_bco].push_back(newhit.release());
  }
}

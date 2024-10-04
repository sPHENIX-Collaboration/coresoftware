#ifndef __TpcTimeFrameBuilder_H__
#define __TpcTimeFrameBuilder_H__

#include <algorithm>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <utility>

class Packet;
class TpcRawHit;
class TH2I;
class TH1;

class TpcTimeFrameBuilder
{
 public:
  explicit TpcTimeFrameBuilder(const int packet_id);
  virtual ~TpcTimeFrameBuilder();

  int ProcessPacket(Packet *);
  bool isMoreDataRequired() const;
  void CleanupUsedPackets(const uint64_t bclk);

  void setVerbosity(const int i) { m_verbosity = i; }

 protected:
  // Length for the 256-bit wide Round Robin Multiplexer for the data stream
  static const size_t DAM_DMA_WORD_LENGTH = 16;

  static const uint16_t FEE_PACKET_MAGIC_KEY_1 = 0xfe;
  static const uint16_t FEE_PACKET_MAGIC_KEY_2 = 0xed;

  static const uint16_t FEE_MAGIC_KEY = 0xba00;
  static const uint16_t GTM_MAGIC_KEY = 0xbb00;
  static const uint16_t GTM_LVL1_ACCEPT_MAGIC_KEY = 0xbbf0;
  static const uint16_t GTM_ENDAT_MAGIC_KEY = 0xbbf1;
  static const uint16_t GTM_MODEBIT_MAGIC_KEY = 0xbbf2;

  static const uint16_t MAX_FEECOUNT = 26;      // that many FEEs
  static const uint16_t MAX_CHANNELS = 8 * 32;  // that many channels per FEE
                                                //  static const uint16_t  HEADER_LENGTH  = 5;
  static const uint16_t HEADER_LENGTH = 7;
  static const uint16_t MAX_PACKET_LENGTH = 1025;

  uint16_t reverseBits(const uint16_t x) const;
  uint16_t crc16(const uint32_t fee, const uint32_t index, const int l) const;
  uint16_t check_data_parity(const unsigned int fee, const unsigned int index, const int l) const;

  //! DMA word structure
  struct dma_word
  {
    uint16_t dma_header;
    uint16_t data[DAM_DMA_WORD_LENGTH - 1];
  };

  int decode_gtm_data(const dma_word &gtm_word);
  int process_fee_data(unsigned int fee_id);

  struct gtm_payload
  {
    uint16_t pkt_type = 0;
    bool is_endat = false;
    bool is_lvl1 = false;
    bool is_modebit = false;
    uint64_t bco = 0;
    uint32_t lvl1_count = 0;
    uint32_t endat_count = 0;
    uint64_t last_bco = 0;
    uint8_t modebits = 0;
    uint8_t userbits = 0;

  };

  struct fee_payload
  {
    uint16_t adc_length = 0;  
    uint16_t data_parity = 0;
    uint16_t sampa_address = 0;
    uint16_t sampa_channel = 0;
    uint16_t channel = 0;
    uint16_t type = 0;
    uint16_t user_word = 0;
    uint16_t bx_timestamp = 0;

    uint16_t data_crc = 0;
    uint16_t calc_crc = 0;
    
    std::vector< std::pair< uint16_t , std::vector<uint16_t> > > waveforms;
  };

  std::vector<std::deque<uint16_t>> m_feeData;

  int m_verbosity = 0;
  int m_packet_id = 0;

  //! common prefix for QA histograms
  std::string m_HistoPrefix;

  //! GTM BCO -> TpcRawHit
  std::map<uint64_t, std::vector<TpcRawHit *>> m_timeFrameMap;
  static const size_t kMaxRawHitLimit = 10000;  // 10k hits per event > 256ch/fee * 26fee

  // -------------------------
  // GTM Matcher
  // -------------------------

  //! GTM Matcher Strategy, order by reliability
  enum enu_gtmMatcherStrategy
  {
    //! if knowing nothing, simply matching FEE waveform to the last level1 tagger
    kLastLv1Tagger = 0,
    //! tracking FEE BCO to GTM BCO sync
    kFEEWaveformBCOSync = 1,
    //! tracking FEE heartbeat to GTM BCO sync
    kFEEHeartBeatSync = 2
  };
  enu_gtmMatcherStrategy m_gtmMatcherStrategy = kLastLv1Tagger;
  uint64_t matchFEE2GTMBCO(uint16_t fee_bco);

  const static int GTMBCObits = 40;
  const static uint64_t GTMBCOmask_ValidBits = (1ULL << GTMBCObits) - 1;
  const static uint64_t GTMBCOmask_RollOverCounts = std::numeric_limits<uint64_t>::max() - GTMBCOmask_ValidBits;
  uint64_t m_GTMBCORollOverCounter = 0;
  uint64_t m_GTMBCOLastReading = 0;
  //! roll over corrected GTM BCO -> GTM payload data
  std::map<uint64_t, gtm_payload> m_gtmData;
  // //! FEE ID -> last roll over corrected GTM BCO
  // std::map<unsigned int, uint64_t> m_feeLastGTMBCO;

  //! errors allowed in match FEE BCO to GTM BCO
  static const int kFEEBCOMatchWindow = 5;
  //! time used to transmit all FEE data to PCIe in GTM BCO time
  static const int kFEEDataTransmissionWindow = 1000000;  // 100ms for very large non-ZS data at 10Hz

  TH2I *m_hFEEDataStream = nullptr;
  TH1 *h_PacketLength = nullptr;
  TH1 *h_PacketLength_Padding = nullptr;
  TH1 *h_PacketLength_Residual = nullptr;
};

#endif

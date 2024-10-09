#ifndef __TpcTimeFrameBuilder_H__
#define __TpcTimeFrameBuilder_H__

#include <algorithm>
#include <cstdint>
#include <deque>
#include <functional>
#include <set>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <list>
#include <optional>

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

  void setVerbosity(const int i);
  
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
    uint32_t bx_timestamp = 0;
    uint64_t gtm_bco = 0;

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
  // Initially developped by Hugo Pereira Da Costa as `MicromegasBcoMatchingInformation`
  // -------------------------
  class BcoMatchingInformation
  {
  public:
    //! constructor
    BcoMatchingInformation() = default;

    //!@name accessor
    //@{

    //! verbosity
    int verbosity() const
    {
      return m_verbosity;
    }

    //! true if matching information is verified
    /**
     * matching information is verified if at least one match
     * between gtm_bco and fee_bco is found
     */
    bool is_verified() const
    {
      return m_verified_from_modebits || m_verified_from_data;
    }

    //! get predicted fee_bco from gtm_bco
    std::optional<uint32_t> get_predicted_fee_bco(uint64_t) const;

    //! multiplier
    double get_gtm_clock_multiplier()
    {
      return m_multiplier;
    }

    //! print gtm bco information
    void print_gtm_bco_information() const;

    //@}

    //!@name modifiers
    //@{

    //! verbosity
    void set_verbosity(int value)
    {
      m_verbosity = value;
    }

    /// set gtm clock multiplier
    void set_gtm_clock_multiplier(double value)
    {
      m_multiplier = value;
    }

    //! find reference from data
    bool find_reference_heartbeat(const fee_payload & HeartBeatPacket);

    //! save all GTM BCO clocks from packet data
    void save_gtm_bco_information(const gtm_payload & gtm_tagger);

    //! find gtm bco matching a given fee
    std::optional<uint64_t> find_gtm_bco(uint32_t /*fee_gtm*/);

    //! cleanup
    void cleanup();

    //! cleanup
    void cleanup(uint64_t /*ref_bco*/);

    //@}

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

    /* see: https://git.racf.bnl.gov/gitea/Instrumentation/sampa_data/src/branch/fmtv2/README.md */
    enum ModeBitType
    {
      BX_COUNTER_SYNC_T = 0,
      ELINK_HEARTBEAT_T = 1,
      SAMPA_EVENT_TRIGGER_T = 2,
      CLEAR_LV1_LAST_T = 6,
      CLEAR_LV1_ENDAT_T = 7
    };

  private:

    //! update multiplier adjustment
    void update_multiplier_adjustment(uint64_t /* gtm_bco */, uint32_t /* fee_bco */);

    //! get adjusted multiplier
    double get_adjusted_multiplier() const;

    //! verbosity
    unsigned int m_verbosity = 0;

    //! verified
    bool m_verified_from_modebits = false;

    bool m_verified_from_data = false;

    //! first lvl1 bco (40 bits)
    uint64_t m_gtm_bco_first = 0;

    //! first fee bco (20 bits)
    uint32_t m_fee_bco_first = 0;

    //! list of available bco
    std::list<uint64_t> m_gtm_bco_list;

    //! matching between fee bco and lvl1 bco
    using m_bco_matching_pair_t = std::pair<unsigned int, uint64_t>;
    std::list<m_bco_matching_pair_t> m_bco_matching_list;

    //! keep track or  fee_bco for which no gtm_bco is found
    std::set<uint32_t> m_orphans;

    //! adjustment to multiplier
    double m_multiplier_adjustment = 0;

    //! running numerator for multiplier adjustment
    double m_multiplier_adjustment_numerator = 0;

    //! running denominator for multiplier adjustment
    double m_multiplier_adjustment_denominator = 0;

    //! running count for multiplier adjustment
    unsigned int m_multiplier_adjustment_count = 0;



    // define limit for matching two fee_bco
    static constexpr unsigned int m_max_multiplier_adjustment_count = 1000;

    // define limit for matching two fee_bco
    static constexpr unsigned int m_max_fee_bco_diff = 10;

    // define limit for matching gtm_bco from lvl1 to enddat

    // define limit for matching fee_bco to fee_bco_predicted
    static constexpr unsigned int m_max_gtm_bco_diff = 100;

  //   // needed to avoid memory leak. Assumes that we will not be assembling more than 50 events at the same time
    static constexpr unsigned int m_max_matching_data_size = 50;

    //! copied from micromegas/MicromegasDefs.h, not available here
    static constexpr int m_nchannels_fee = 256;

    // get the difference between two BCO.
    template <class T>
    inline static constexpr T get_bco_diff(const T& first, const T& second)
    {
      return first < second ? (second - first) : (first - second);
    }

    // this is the clock multiplier from lvl1 to fee clock
    // Tested with Run24 data. Could be changable in future runs
    double m_multiplier = 4.262916255;

  };

  //! map bco_information_t to packet id
  std::vector<BcoMatchingInformation> m_bcoMatchingInformation_vec;

  TH2I *m_hFEEDataStream = nullptr;
  TH1 *h_PacketLength = nullptr;
  TH1 *h_PacketLength_Padding = nullptr;
  TH1 *h_PacketLength_Residual = nullptr;
};

#endif

#ifndef Fun4All_TpcTimeFrameBuilder_H
#define Fun4All_TpcTimeFrameBuilder_H

#include <algorithm>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

class Packet;
class TpcRawHit;
class PHTimer;
class TH1;
class TH2;
class TTree;

// NOLINTNEXTLINE(hicpp-special-member-functions)
class TpcTimeFrameBuilder
{
 public:
  explicit TpcTimeFrameBuilder(const int packet_id);
  virtual ~TpcTimeFrameBuilder();

  int ProcessPacket(Packet *);
  bool isMoreDataRequired(const uint64_t &gtm_bco) const;
  void CleanupUsedPackets(const uint64_t &bclk);
  std::vector<TpcRawHit *> &getTimeFrame(const uint64_t &gtm_bco);

  void setVerbosity(const int i);
  void setFastBCOSkip(bool fastBCOSkip = true)
  {
    m_fastBCOSkip = fastBCOSkip;
  }

  // enable saving of digital current debug TTree with file name `name`
  void SaveDigitalCurrentDebugTTree(const std::string &name);

 protected:
  // Length for the 256-bit wide Round Robin Multiplexer for the data stream
  static const size_t DAM_DMA_WORD_LENGTH = 16;

  static const uint16_t FEE_PACKET_MAGIC_KEY_1 = 0xfe;
  static const uint16_t FEE_PACKET_MAGIC_KEY_2 = 0xed;
  static const uint16_t FEE_PACKET_MAGIC_KEY_3_DC = 0xdcdc; // Digital Current word[3]

  static const uint16_t FEE_MAGIC_KEY = 0xba00;
  static const uint16_t GTM_MAGIC_KEY = 0xbb00;
  static const uint16_t GTM_LVL1_ACCEPT_MAGIC_KEY = 0xbbf0;
  static const uint16_t GTM_ENDAT_MAGIC_KEY = 0xbbf1;
  static const uint16_t GTM_MODEBIT_MAGIC_KEY = 0xbbf2;

  static const uint16_t MAX_FEECOUNT = 26;              // that many FEEs
  static const uint16_t MAX_SAMPA = 8;                  // that many FEEs
  static const uint16_t MAX_CHANNELS = MAX_SAMPA * 32;  // that many channels per FEE
                                                        //  static const uint16_t  HEADER_LENGTH  = 5;
  static const uint16_t HEADER_LENGTH = 7;
  static const uint16_t MAX_PACKET_LENGTH = 1025;

  static const uint16_t GL1_BCO_MATCH_WINDOW = 256;  // BCOs

  int m_hitFormat = -1;

  uint16_t reverseBits(const uint16_t x) const;
  std::pair<uint16_t, uint16_t> crc16_parity(const uint32_t fee, const uint16_t l) const;

  //! DMA word structure
  struct dma_word
  {
    uint16_t dma_header = 0;
    uint16_t data[DAM_DMA_WORD_LENGTH - 1] = {0};
  };

  int decode_gtm_data(const dma_word &gtm_word);
  int process_fee_data(unsigned int fee_id);
  void process_fee_data_waveform(const unsigned int & fee_id, std::deque<uint16_t>& data_buffer);
  void process_fee_data_digital_current(const unsigned int & fee_id, std::deque<uint16_t>& data_buffer);

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
    uint16_t fee_id = 0;
    uint16_t adc_length = 0;
    uint16_t sampa_address = 0;
    uint16_t sampa_channel = 0;
    uint16_t channel = 0;
    uint16_t type = 0;
    uint16_t user_word = 0;
    uint32_t bx_timestamp = 0;
    uint64_t gtm_bco = 0;

    uint16_t data_crc = 0;
    uint16_t calc_crc = 0;
    
    uint16_t data_parity = 0;
    uint16_t calc_parity = 0;

    std::vector<std::pair<uint16_t, std::vector<uint16_t>>> waveforms;
  };

  struct digital_current_payload
  {
    static const int MAX_CHANNELS = 8;

    uint16_t fee {std::numeric_limits<uint16_t>::max()};
    uint16_t pkt_length {std::numeric_limits<uint16_t>::max()};
    uint16_t channel {std::numeric_limits<uint16_t>::max()};
    // uint16_t sampa_max_channel {std::numeric_limits<uint16_t>::max()};
    uint16_t sampa_address {std::numeric_limits<uint16_t>::max()};
    uint32_t bx_timestamp {0};
    uint32_t current[MAX_CHANNELS] {0};
    uint32_t nsamples[MAX_CHANNELS] {0};
    uint16_t data_crc {std::numeric_limits<uint16_t>::max()};
    uint16_t calc_crc = {std::numeric_limits<uint16_t>::max()};
    // uint16_t type {std::numeric_limits<uint16_t>::max()};
  };

  class DigitalCurrentDebugTTree
  {
   public:
    explicit DigitalCurrentDebugTTree(const std::string &name);
    virtual ~DigitalCurrentDebugTTree();

    void fill(const digital_current_payload &payload);

   private:
    digital_current_payload m_payload;

    std::string m_name;
    TTree *m_tDigitalCurrent = nullptr;
  };
  DigitalCurrentDebugTTree * m_digitalCurrentDebugTTree = nullptr;

  // -------------------------
  // GTM Matcher
  // Initially developped by Hugo Pereira Da Costa as `MicromegasBcoMatchingInformation`
  // -------------------------
  class BcoMatchingInformation
  {
   public:
    //! constructor
    explicit BcoMatchingInformation(const std::string &name);

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

    //! matching between fee bco and lvl1 bco
    using m_gtm_fee_bco_matching_pair_t = std::pair<uint64_t, uint32_t>;
    using m_fee_gtm_bco_matching_pair_t = std::pair<uint32_t, uint64_t>;

    //! get reference bco
    const std::optional<m_gtm_fee_bco_matching_pair_t> &get_reference_bco() const
    {
      return m_bco_reference;
    }

    //! whether FEE data has moved pass the given gtm_bco
    bool isMoreDataRequired(const uint64_t &gtm_bco) const;

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

    /// set gtm clock with rollover correction
    uint64_t get_gtm_rollover_correction(const uint64_t &gtm_bco) const;

    //! find reference from data
    std::optional<uint64_t> find_reference_heartbeat(const fee_payload &HeartBeatPacket);

    //! save all GTM BCO clocks from packet data
    void save_gtm_bco_information(const gtm_payload &gtm_tagger);

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

    // Command   | OLD Mode-Bit | New Mode-Number | Function
    // ======================================================
    // NOP       |     0b000    |             0x0 | No Operation
    // BX_SYNC   |     0b001    |             0x1 | SAMPA Beam-crossing sync
    // H_BEAT    |     0b010    |             0x2 | Generates Heartbeat frame
    // TRIG      |     0b100    |             0x3 | Trigger data when FEM user bit is 0b01, otherwise the level 1 accept is used when FEM user bit is 0b00
    // CLK_SYNC  |     N/A      |             0x4 | Reset and align 40 MHz and 20 MHz clocks to SAMPA
    // SAMPA_RST |     N/A      |             0x5 | Hard reset SAMPA
    // DC_START  |     N/A      |             0x6 | Start digital current reading
    // DC_STOP   |     N/A      |             0x7 | Stop and send digital current packet
    enum ModeBitType
    {
      BX_COUNTER_SYNC_T = 0x1,
      ELINK_HEARTBEAT_T = 0x2
      // SAMPA_EVENT_TRIGGER_T = 2,
      // CLEAR_LV1_LAST_T = 6,
      // CLEAR_LV1_ENDAT_T = 7
    };

    // get the difference between two BCO WITHOUT rollover corrections
    template <class T>
    inline static constexpr T get_bco_diff(
        const T &first, const T &second)
    {
      return first < second ? (second - first) : (first - second);
    }

    // get the difference between two BCO with rollover corrections
    inline static constexpr uint32_t get_fee_bco_diff(
        const uint32_t &first, const uint32_t &second)  // NOLINT(misc-unused-parameters)
    {
      const uint32_t diff_raw = get_bco_diff(first, second);

      return (diff_raw < (1U << (m_FEE_CLOCK_BITS / 2))) ? diff_raw : (1U << m_FEE_CLOCK_BITS) - diff_raw;
    }

   private:
    std::string m_name;

    //! verbosity
    unsigned int m_verbosity = 0;

    //! verified
    bool m_verified_from_modebits = false;

    bool m_verified_from_data = false;

    //! list of available bco, sorted in time with rollover corrected
    std::list<uint64_t> m_gtm_bco_trig_list;

    //! list of available GTM -> FEE bco mapping for synchronization
    std::optional<m_gtm_fee_bco_matching_pair_t> m_bco_reference = std::nullopt;

    // std::optional< std::pair< uint64_t, uint32_t > > m_bco_reference_candidate = std::nullopt;
    //! not yet matched heart beats
    std::list<m_gtm_fee_bco_matching_pair_t> m_bco_reference_candidate_list;
    static constexpr unsigned int m_max_bco_reference_candidate_list_size = 16;

    // //! list of heart beat GTM BCO that is still to be matched
    // std::queue<uint64_t> m_heartbeat_gtm_bco_queue;
    // static constexpr unsigned int m_max_heartbeat_queue_size = 16;

    //! list of available GTM -> FEE bco mapping for trigger association
    std::map<uint64_t, uint32_t> m_gtm_bco_trigger_map;

    std::list<m_fee_gtm_bco_matching_pair_t> m_bco_matching_list;

    //! keep track or  fee_bco for which no gtm_bco is found
    std::set<uint32_t> m_orphans;

    // define limit for matching two lvl1 and EnDAT tagger BCOs
    static constexpr int m_max_lv1_endat_bco_diff = 16;

    // define limit for matching two fee_bco
    static constexpr unsigned int m_max_fee_bco_diff = 64;

    // define limit for matching gtm_bco from lvl1 to enddat

    // define limit for matching fee_bco to fee_bco_predicted
    static constexpr unsigned int m_max_gtm_bco_diff = 256;

    //   // needed to avoid memory leak. Assumes that we will not be assembling more than 50 events at the same time
    static constexpr unsigned int m_max_matching_data_size = 10;

    //! max time in GTM BCO for FEE data to sync over to datastream
    static constexpr unsigned int m_max_fee_sync_time = 1024 * 8;

    static constexpr unsigned int m_FEE_CLOCK_BITS = 20;
    static constexpr unsigned int m_GTM_CLOCK_BITS = 40;

    double m_multiplier = 0;

    TH1 *m_hNorm = nullptr;
    TH1 *m_hFEEClockAdjustment_MatchedReference = nullptr;
    TH1 *m_hFEEClockAdjustment_MatchedNew = nullptr;
    TH1 *m_hFEEClockAdjustment_Unmatched = nullptr;
    TH1 *m_hGTMNewEventSpacing = nullptr;
    TH1 *m_hFindGTMBCO_MatchedExisting_BCODiff = nullptr;
    TH1 *m_hFindGTMBCO_MatchedNew_BCODiff = nullptr;

  };  //   class BcoMatchingInformation

 private:
  std::vector<std::deque<uint16_t>> m_feeData;

  int m_verbosity = 0;
  int m_packet_id = 0;

  //! common prefix for QA histograms
  std::string m_HistoPrefix;

  //! GTM BCO -> TpcRawHit
  //! Map to store TpcRawHit pointers indexed by GTM BCO values
  //! This is used to organize hits into time frames based on their BCO values
  std::map<uint64_t, std::vector<TpcRawHit *>> m_timeFrameMap;
  static const size_t kMaxRawHitLimit = 10000;  // 10k hits per event > 256ch/fee * 26fee
  std::queue<uint64_t> m_UsedTimeFrameSet;

  //! fast skip mode when searching for particular GL1 BCO over long segment of files
  bool m_fastBCOSkip = false;

  //! map bco_information_t to packet id
  std::vector<BcoMatchingInformation> m_bcoMatchingInformation_vec;

  //! QA area

  PHTimer *m_packetTimer = nullptr;

  TH1 *m_hNorm = nullptr;
  TH2 *m_hFEEDataStream = nullptr;
  TH1 *m_hFEEChannelPacketCount = nullptr;
  TH2 *m_hFEESAMPAADC = nullptr;
  TH1 *m_hFEESAMPAHeartBeatSync = nullptr;

  TH1 *h_PacketLength = nullptr;
  TH1 *h_PacketLength_Padding = nullptr;
  TH1 *h_PacketLength_Residual = nullptr;

  TH1 *h_GTMClockDiff_Matched = nullptr;
  TH1 *h_GTMClockDiff_Unmatched = nullptr;
  TH1 *h_GTMClockDiff_Dropped = nullptr;
  TH1 *h_TimeFrame_Matched_Size = nullptr;

  TH2 *h_ProcessPacket_Time = nullptr;
};

#endif

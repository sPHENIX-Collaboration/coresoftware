#ifndef FUN4ALLRAW_MICROMEGASBCOMATCHINGINFORMATION_V2_H
#define FUN4ALLRAW_MICROMEGASBCOMATCHINGINFORMATION_V2_H

/*!
 * \file MicromegasBcoMatchingInformation_v2.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief handles matching between GTM and FEE BCO clock
 */

#include <cstdint>
#include <list>
#include <optional>
#include <set>
#include <vector>
#include <utility>

class MicromegasBcoMatchingInformation_v2
{
 public:

  //! constructor
  MicromegasBcoMatchingInformation_v2() = default;

  //! gtm data
  class gtm_payload
  {
    public:

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

  //! fee data
  class fee_payload
  {
    public:
    uint16_t adc_length = 0;
    uint16_t data_parity = 0;
    uint16_t sampa_address = 0;
    uint16_t sampa_channel = 0;
    uint16_t channel = 0;
    uint16_t type = 0;
    uint16_t user_word = 0;
    uint32_t bx_timestamp = 0;
    uint16_t data_crc = 0;
    uint16_t calc_crc = 0;

    std::vector< std::pair< uint16_t , std::vector<uint16_t> > > waveforms;
  };

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
  static bool gtm_clock_multiplier_is_set()
  {
    return m_multiplier_is_set;
  }

  //! multiplier
  static double get_gtm_clock_multiplier()
  {
    return m_multiplier;
  }

  //! get adjusted multiplier
  double get_adjusted_multiplier() const;

  //! print gtm bco information
  void print_gtm_bco_information() const;

  //! get BCO matching reference
  using m_bco_matching_pair_t = std::pair<uint32_t, uint64_t>;
  const m_bco_matching_pair_t& get_bco_matching_reference() const
  { return m_bco_reference; }

  //! get last gtm bco
  uint64_t get_gtm_bco_last() const
  { return m_gtm_bco_list.empty() ? 0:*m_gtm_bco_list.rbegin(); }

  //! returns true if more data needs to be fetched.
  /** it is based on the latest BCO read from the data stream, from either tagger or heartbeat */
  bool is_more_data_required(uint64_t /*gtm_bco*/) const;

  //@}

  //!@name modifiers
  //@{

  //! verbosity
  void set_verbosity(int value)
  {
    m_verbosity = value;
  }

  /// set gtm clock multiplier
  static void set_gtm_clock_multiplier(double value)
  {
    m_multiplier_is_set = true;
    m_multiplier = value;
  }

  /// enable multiplier adjustment
  static void set_enable_multiplier_adjustment( bool value )
  { m_multiplier_adjustment_enabled = value; }

  /// muliplier adjustment count
  /** controls how often the gtm multiplier is automatically adjusted */
  static void set_max_multiplier_adjustment_count( unsigned int value )
  {
    m_max_multiplier_adjustment_count = value;
  }

  // define limit for matching fee_bco to fee_bco_predicted from gtm_bco
  static void set_gtm_bco_diff( unsigned int value )
  {
    m_max_gtm_bco_diff = value;
  }

  //! max time in GTM BCO for FEE data to sync over to datastream
  static void set_m_max_fee_sync_time( unsigned int value )
  {
    m_max_fee_sync_time = value;
  }

  //! find reference from modebits
  bool find_reference_from_modebits(const gtm_payload&);

  //! find reference from data
  bool find_reference_from_data(const fee_payload&);

  //! save all GTM BCO clocks from packet data
  void save_gtm_bco_information(int /* packet_id */, const gtm_payload&);

  //! find gtm bco matching a given fee
  /**
   * packet and fee ids are not necessary to the calculation.
   * They are pased here for the clarity of debugging messages
   */
  std::optional<uint64_t> find_gtm_bco(int /*packet_id*/, unsigned int /*fee_id*/, uint32_t /*fee_gtm*/);

  //! cleanup
  void cleanup();

  //! cleanup
  void cleanup(uint64_t /*ref_bco*/);

  //@}

  //!@name utilities
  //@{

  //! get difference between two GTM BCO, properly accounting for 40bits rollover
  /** based on Jin's code in TpcTimeFrameBuilder */
  static int64_t get_signed_gtm_bco_diff(uint64_t /*first*/, uint64_t /*second*/);

  //! get difference between two FEE BCO, properly accounting for 20bits rollover
  /** based on Jin's code in TpcTimeFrameBuilder */
  static int32_t get_signed_fee_bco_diff(uint32_t /*first*/, uint32_t /*second*/);

  //! get difference between two GTM BCO, properly accounting for 40bits rollover
  static uint64_t get_unsigned_gtm_bco_diff(uint64_t /*first*/, uint64_t /*second*/);

  //! get difference between two FEE BCO, properly accounting for 20bits rollover
  static uint32_t get_unsigned_fee_bco_diff(uint32_t /*first*/, uint32_t /*second*/);

  //@}

 private:

  //! update multiplier adjustment
  void update_multiplier_adjustment(uint64_t /* gtm_bco */, uint32_t /* fee_bco */);

  //! verbosity
  unsigned int m_verbosity = 0;

  //! verified
  bool m_verified_from_modebits = false;
  bool m_verified_from_data = false;

  //! last found fee_bco
  /** used to try finding bco reference from data */
  uint32_t m_fee_bco_prev = 0;
  bool m_has_fee_bco_prev = false;

  //! list of available bco
  std::list<uint64_t> m_gtm_bco_list;

  //! reference matching
  m_bco_matching_pair_t m_bco_reference;

  //! matching between fee bco and lvl1 bco
  std::list<m_bco_matching_pair_t> m_bco_matching_list;

  //! keep track or  fee_bco for which no gtm_bco is found
  std::set<uint32_t> m_orphans;

  //! true if multiplier is set
  static bool m_multiplier_is_set;

  //! gtm clock multiplier
  static double m_multiplier;

  //! true if multiplier adjustment is enabled
  static bool m_multiplier_adjustment_enabled;

  //! multiplier adjustment count
  /* controls how often the gtm multiplier is automatically adjusted */
  static unsigned int m_max_multiplier_adjustment_count;

  // define limit for matching fee_bco to fee_bco_predicted
  static unsigned int m_max_gtm_bco_diff;

  //! max time in GTM BCO for FEE data to sync over to datastream
  static unsigned int m_max_fee_sync_time;

  //! adjustment to multiplier
  double m_multiplier_adjustment = 0;

  //! running numerator for multiplier adjustment
  double m_multiplier_adjustment_numerator = 0;

  //! running denominator for multiplier adjustment
  double m_multiplier_adjustment_denominator = 0;

  //! running count for multiplier adjustment
  unsigned int m_multiplier_adjustment_count = 0;
};

#endif

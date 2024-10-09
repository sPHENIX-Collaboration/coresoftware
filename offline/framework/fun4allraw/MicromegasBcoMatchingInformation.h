#ifndef MICROMEGAS_MicromegasBcoMatchingInformation_H
#define MICROMEGAS_MicromegasBcoMatchingInformation_H

/*!
 * \file MicromegasBcoMatchingInformation.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief handles matching between GTM and FEE BCO clock
 */

#include <cstdint>
#include <list>
#include <optional>
#include <set>
#include <utility>

class Packet;

class MicromegasBcoMatchingInformation
{
 public:
  //! constructor
  MicromegasBcoMatchingInformation() = default;

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
  static double get_gtm_clock_multiplier()
  {
    return m_multiplier;
  }

  //! get adjusted multiplier
  double get_adjusted_multiplier() const;

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
  static void set_gtm_clock_multiplier(double value)
  {
    m_multiplier = value;
  }

  //! find reference from modebits
  bool find_reference_from_modebits(Packet*);

  //! find reference from data
  bool find_reference_from_data(Packet*);

  //! save all GTM BCO clocks from packet data
  void save_gtm_bco_information(Packet*);

  //! find gtm bco matching a given fee
  std::optional<uint64_t> find_gtm_bco(uint32_t /*fee_gtm*/);

  //! cleanup
  void cleanup();

  //! cleanup
  void cleanup(uint64_t /*ref_bco*/);

  //@}

 private:

  //! update multiplier adjustment
  void update_multiplier_adjustment(uint64_t /* gtm_bco */, uint32_t /* fee_bco */);

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

  //! gtm clock multiplier
  static double m_multiplier;

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

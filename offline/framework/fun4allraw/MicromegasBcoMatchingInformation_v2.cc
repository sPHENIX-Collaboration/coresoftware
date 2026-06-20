
/*!
 * \file MicromegasBcoMatchingInformation_v2.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief handles matching between GTM and FEE BCO clock
 */

#include "MicromegasBcoMatchingInformation_v2.h"

#include <phool/phool.h>

#include <Event/packet.h>

#include <TSystem.h>

#include <algorithm>
#include <vector>

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

  // define limit for matching two fee_bco
  constexpr uint32_t m_max_fee_bco_diff = 10;

  // needed to avoid memory leak. Assumes that we will not be assembling more than 50 events at the same time
  constexpr unsigned int m_max_matching_data_size = 50;

  //! copied from micromegas/MicromegasDefs.h, not available here
  constexpr int m_nchannels_fee = 256;

  // gtm clock bits
  /* used for rollover calculation */
  constexpr unsigned int m_GTM_CLOCK_BITS = 40U;
  constexpr uint64_t m_GTM_CLOCK_MASK = (1ULL << m_GTM_CLOCK_BITS)-1;
  constexpr int64_t m_GTM_CLOCK_RANGE = 1ULL << m_GTM_CLOCK_BITS;
  constexpr int64_t m_GTM_CLOCK_HALF_RANGE = 1ULL << (m_GTM_CLOCK_BITS-1);

  // Fee clock bits
  constexpr unsigned int m_FEE_CLOCK_BITS = 20U;
  constexpr uint32_t m_FEE_CLOCK_MASK = (1UL << m_FEE_CLOCK_BITS)-1ULL;
  constexpr int32_t m_FEE_CLOCK_RANGE = 1UL << m_FEE_CLOCK_BITS;
  constexpr int32_t m_FEE_CLOCK_HALF_RANGE = 1UL << (m_FEE_CLOCK_BITS-1UL);

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
    BX_COUNTER_SYNC_T = 0b001,
    ELINK_HEARTBEAT_T = 0b010
  };

}  // namespace

//! this is the clock multiplier from lvl1 to fee clock
bool MicromegasBcoMatchingInformation_v2::m_multiplier_is_set = false;
double MicromegasBcoMatchingInformation_v2::m_multiplier = 0;

//! true if on-fly multiplier adjustment is enabled
bool MicromegasBcoMatchingInformation_v2::m_multiplier_adjustment_enabled = true;

//! muliplier adjustment count
/* controls how often the gtm multiplier is automatically adjusted */
unsigned int MicromegasBcoMatchingInformation_v2::m_max_multiplier_adjustment_count = 200;

//! define limit for matching fee_bco to fee_bco_predicted
unsigned int MicromegasBcoMatchingInformation_v2::m_max_gtm_bco_diff = 60;

//! Max time forward to ensure that a given
unsigned int MicromegasBcoMatchingInformation_v2::m_max_fee_sync_time = 1024 * 96;

//___________________________________________________
int64_t MicromegasBcoMatchingInformation_v2::get_signed_gtm_bco_diff(uint64_t first, uint64_t second)
{
  // calculate raw diff
  int64_t diff = static_cast<int64_t>(first & m_GTM_CLOCK_MASK) - static_cast<int64_t>(second & m_GTM_CLOCK_MASK);

  // make sure result is within +/- m_GTM_CLOCK_HALF_RANGE
  if (diff > m_GTM_CLOCK_HALF_RANGE) { diff -= m_GTM_CLOCK_RANGE; }
  else if (diff < -m_GTM_CLOCK_HALF_RANGE) { diff += m_GTM_CLOCK_RANGE; }

  return diff;
}

//___________________________________________________
int32_t MicromegasBcoMatchingInformation_v2::get_signed_fee_bco_diff(uint32_t first, uint32_t second)
{
  // calculate raw diff
  int32_t diff = static_cast<int32_t>(first & m_FEE_CLOCK_MASK) - static_cast<int32_t>(second & m_FEE_CLOCK_MASK);

  // make sure result is within +/- m_FEE_CLOCK_HALF_RANGE
  if (diff > m_FEE_CLOCK_HALF_RANGE) { diff -= m_FEE_CLOCK_RANGE; }
  else if (diff < -m_FEE_CLOCK_HALF_RANGE) { diff += m_FEE_CLOCK_RANGE; }

  return diff;
}

//___________________________________________________
uint64_t MicromegasBcoMatchingInformation_v2::get_unsigned_gtm_bco_diff(uint64_t first, uint64_t second)
{
  const auto diff = get_signed_gtm_bco_diff(first, second);
  return uint64_t((diff<0) ? -diff:diff);
}

//___________________________________________________
uint32_t MicromegasBcoMatchingInformation_v2::get_unsigned_fee_bco_diff(uint32_t first, uint32_t second)
{
  const auto diff = get_signed_fee_bco_diff(first, second);
  return uint32_t((diff<0) ? -diff:diff);
}

//___________________________________________________
std::optional<uint32_t> MicromegasBcoMatchingInformation_v2::get_predicted_fee_bco(uint64_t gtm_bco) const
{
  // check proper initialization
  if (!is_verified())
  {
    return std::nullopt;
  }

  // get gtm bco difference
  const int64_t gtm_bco_difference = get_signed_gtm_bco_diff(gtm_bco, m_bco_reference.second);

  // convert to fee bco, and truncate to 20 bits
  const int64_t fee_bco_predicted = m_bco_reference.first + get_adjusted_multiplier() * gtm_bco_difference;
  return static_cast<uint32_t>(fee_bco_predicted) & m_FEE_CLOCK_MASK;
}

//___________________________________________________
void MicromegasBcoMatchingInformation_v2::print_gtm_bco_information() const
{
  if (!m_gtm_bco_list.empty())
  {
    std::cout
        << "MicromegasBcoMatchingInformation_v2::print_gtm_bco_information -"
        << " gtm_bco: " << std::hex << m_gtm_bco_list << std::dec
        << std::endl;

    // also print predicted fee bco
    if (is_verified())
    {
      std::list<uint32_t> fee_bco_predicted_list;
      std::transform(
          m_gtm_bco_list.begin(),
          m_gtm_bco_list.end(),
          std::back_inserter(fee_bco_predicted_list),
          [this](const uint64_t& gtm_bco)
          { return get_predicted_fee_bco(gtm_bco).value(); });

      std::cout
          << "MicromegasBcoMatchingInformation_v2::print_gtm_bco_information -"
          << " fee_bco_predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
    }
  }
}

//___________________________________________________
bool MicromegasBcoMatchingInformation_v2::is_more_data_required( uint64_t gtm_bco ) const
{
  // check proper initialization
  if( !is_verified() ) { return true; }

  // compare to reference
  if( get_signed_gtm_bco_diff( m_bco_reference.second, gtm_bco ) > m_max_fee_sync_time )
  { return false; }

  // check against stored bco
  if( !m_gtm_bco_list.empty() )
  {
    if( get_signed_gtm_bco_diff( m_gtm_bco_list.back(), gtm_bco ) > m_max_fee_sync_time )
    { return false; }
  }

  // check against matched BCOs
  if( !m_bco_matching_list.empty() )
  {
    if( get_signed_gtm_bco_diff( m_bco_matching_list.back().second, gtm_bco ) > m_max_fee_sync_time )
    { return false; }
  }

  return true;
}


//___________________________________________________
void MicromegasBcoMatchingInformation_v2::save_gtm_bco_information(int /*packet_id*/, const MicromegasBcoMatchingInformation_v2::gtm_payload& payload)
{
  if (payload.is_lvl1)
  {
    // save lvl1 BCO
    const auto& gtm_bco = payload.bco;
    m_gtm_bco_list.push_back(gtm_bco);
  }
  else if (payload.is_endat)
  {
    // also save ENDDAT bco
    const auto& gtm_bco = payload.bco;

    // add to list if difference to last entry is big enough
    if (m_gtm_bco_list.empty() || get_signed_gtm_bco_diff(gtm_bco,m_gtm_bco_list.back()) > 10)
    {
      m_gtm_bco_list.push_back(gtm_bco);
    }
  }
  else if (payload.is_modebit)
  {
    // also save hearbeats BCO
    const auto& modebits = payload.modebits;
    if (modebits == ELINK_HEARTBEAT_T)
    {
      const auto& gtm_bco = payload.bco;
      m_gtm_bco_list.push_back(gtm_bco);
    }
  }
}

//___________________________________________________
bool MicromegasBcoMatchingInformation_v2::find_reference_from_modebits(const MicromegasBcoMatchingInformation_v2::gtm_payload& payload)
{
  if (payload.is_modebit)
  {
    // get modebits
    const auto& modebits = payload.modebits;
    if (modebits == BX_COUNTER_SYNC_T)
    {
      std::cout << "MicromegasBcoMatchingInformation_v2::find_reference_from_modebits"
                << " found reference from modebits"
                << std::endl;

      // get BCO and assign
      const auto& gtm_bco = payload.bco;
      m_bco_reference = {0, gtm_bco};
      m_verified_from_modebits = true;
      return true;
    }
  }
  return false;
}

//___________________________________________________
bool MicromegasBcoMatchingInformation_v2::find_reference_from_data(const fee_payload& payload)
{
  // store gtm bco and diff to previous in an array
  std::vector<uint64_t> gtm_bco_list;
  std::vector<int64_t> fee_bco_diff_list;
  for (const auto& gtm_bco : m_gtm_bco_list)
  {
    if (!gtm_bco_list.empty())
    {
      // add difference to last
      // get gtm bco difference with proper rollover accounting
      const int64_t gtm_bco_difference = get_signed_gtm_bco_diff( gtm_bco, gtm_bco_list.back() );
      fee_bco_diff_list.push_back(get_adjusted_multiplier() * gtm_bco_difference);
    }

    // append to list
    gtm_bco_list.push_back(gtm_bco);
  }

  // print all differences
  if (verbosity())
  {
    std::cout << "MicromegasBcoMatchingInformation_v2::find_reference_from_data - fee_bco_diff_list: " << fee_bco_diff_list << std::endl;
  }

  // skip hearbeat
  if (payload.type == HEARTBEAT_T)
  {
    return false;
  }

  // bound check
  if (payload.channel >= m_nchannels_fee)
  {
    return false;
  }

  // get timestamp
  const auto& fee_bco = payload.bx_timestamp;
  if (!m_has_fee_bco_prev)
  {
    m_fee_bco_prev = fee_bco;
    m_has_fee_bco_prev = true;
  }

  // calculate difference
  const uint32_t fee_bco_diff = get_unsigned_fee_bco_diff(fee_bco, m_fee_bco_prev);

  // discard identical fee_bco
  if (fee_bco_diff < m_max_fee_bco_diff)
  {
    return false;
  }

  std::cout << "MicromegasBcoMatchingInformation_v2::find_reference_from_data - fee_bco_diff: " << fee_bco_diff << std::endl;

  // look for matching diff in gtm_bco array
  for (size_t i = 0; i < fee_bco_diff_list.size(); ++i)
  {
    uint32_t sum = 0;
    for (size_t j = i; j < fee_bco_diff_list.size(); ++j)
    {
      sum += fee_bco_diff_list[j];
      if (get_unsigned_fee_bco_diff(sum, fee_bco_diff) < m_max_fee_bco_diff)
      {
        m_verified_from_data = true;
        m_bco_reference = { m_fee_bco_prev, gtm_bco_list[i] };

        if (verbosity())
        {
          std::cout << "MicromegasBcoMatchingInformation_v2::find_reference_from_data - matching is verified" << std::endl;
          std::cout
              << "MicromegasBcoMatchingInformation_v2::find_reference_from_data -"
              << std::hex
              << " m_bco_reference: ( 0x"
              << m_bco_reference.first
              << ", 0x" << m_bco_reference.second << ")"
              << std::dec
              << std::endl;
        }
        return true;
      }
    }
  }

  // update previous fee_bco
  m_fee_bco_prev = fee_bco;
  return false;
}

//___________________________________________________
std::optional<uint64_t> MicromegasBcoMatchingInformation_v2::find_gtm_bco(int packet_id, unsigned int fee_id, uint32_t fee_bco)
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
      [fee_bco](const m_bco_matching_pair_t& pair)
      { return get_unsigned_fee_bco_diff(pair.first, fee_bco) < m_max_fee_bco_diff; });

  if (bco_matching_iter != m_bco_matching_list.end())
  {
    return bco_matching_iter->second;
  }

  // find element for which predicted fee_bco matches fee_bco, within limit
  const auto iter = std::find_if(
      m_gtm_bco_list.begin(),
      m_gtm_bco_list.end(),
      [this, fee_bco](const uint64_t& gtm_bco)
      { return get_unsigned_fee_bco_diff(get_predicted_fee_bco(gtm_bco).value(), fee_bco) < m_max_gtm_bco_diff; });

  // check
  if (iter != m_gtm_bco_list.end())
  {
    const auto gtm_bco = *iter;
    if (verbosity())
    {
      if (auto opt_fee_bco = get_predicted_fee_bco(gtm_bco))  // check if optional exists
      {
        const auto fee_bco_predicted = *opt_fee_bco;  // get_predicted_fee_bco(gtm_bco).value();
        const auto fee_bco_diff = get_unsigned_fee_bco_diff(fee_bco_predicted, fee_bco);

        std::cout << "MicromegasBcoMatchingInformation_v2::find_gtm_bco -"
          << " packet_id: " << packet_id
          << " fee_id: " << fee_id
          << std::hex
          << " fee_bco: 0x" << fee_bco
          << " predicted: 0x" << fee_bco_predicted
          << " gtm_bco: 0x" << gtm_bco
          << std::dec
          << " difference: " << fee_bco_diff
          << std::endl;
      }
    }
    // save fee_bco and gtm_bco matching in map
    m_bco_matching_list.emplace_back(fee_bco, gtm_bco);

    // remove gtm bco from runing list
    m_gtm_bco_list.erase(iter);

    // update clock adjustment
    if (m_multiplier_adjustment_enabled)
    {
      update_multiplier_adjustment(gtm_bco, fee_bco);
    }

    return gtm_bco;
  }

  if (m_orphans.insert(fee_bco).second)
  {
    if (verbosity())
    {
      // find element for which predicted fee_bco is the closest to request
      const auto iter2 = std::min_element(
          m_gtm_bco_list.begin(),
          m_gtm_bco_list.end(),
          [this, fee_bco](const uint64_t& first, const uint64_t& second)
          { return get_unsigned_fee_bco_diff(get_predicted_fee_bco(first).value(), fee_bco) < get_unsigned_fee_bco_diff(get_predicted_fee_bco(second).value(), fee_bco); });

      // compared to the previous statement, this checks if the optional
      int fee_bco_diff = -1;

      if (iter2 != m_gtm_bco_list.end())
      {
        auto predicted = get_predicted_fee_bco(*iter2);

        if (predicted)
        {
          fee_bco_diff = get_unsigned_fee_bco_diff(*predicted, fee_bco);
        }
      }

      std::cout << "MicromegasBcoMatchingInformation_v2::find_gtm_bco -"
                << " packet_id: " << packet_id
                << " fee_id: " << fee_id
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

//___________________________________________________
void MicromegasBcoMatchingInformation_v2::cleanup()
{
  // remove old gtm_bco and matching
  while (m_gtm_bco_list.size() > m_max_matching_data_size)
  {
    m_gtm_bco_list.pop_front();
  }
  while (m_bco_matching_list.size() > m_max_matching_data_size)
  {
    m_bco_matching_list.pop_front();
  }

  // clear orphans
  m_orphans.clear();
}

//___________________________________________________
void MicromegasBcoMatchingInformation_v2::cleanup(uint64_t ref_bco)
{
  // erase all elements from bco_list that are less than or equal to ref_bco
  m_gtm_bco_list.erase(std::remove_if(m_gtm_bco_list.begin(), m_gtm_bco_list.end(),
    [ref_bco](const uint64_t& bco) { return get_signed_gtm_bco_diff( bco,ref_bco ) <= 0; }), m_gtm_bco_list.end());

  // erase all elements from bco_list that are less than or equal to ref_bco
  m_bco_matching_list.erase(std::remove_if(m_bco_matching_list.begin(), m_bco_matching_list.end(),
    [ref_bco](const m_bco_matching_pair_t& pair) { return get_signed_gtm_bco_diff( pair.second, ref_bco ) <= 0; }),
    m_bco_matching_list.end());

  // clear orphans
  m_orphans.clear();
}

//___________________________________________________
double MicromegasBcoMatchingInformation_v2::get_adjusted_multiplier() const
{
  return m_multiplier + m_multiplier_adjustment;
}

//___________________________________________________
void MicromegasBcoMatchingInformation_v2::update_multiplier_adjustment(uint64_t gtm_bco, uint32_t fee_bco)
{
  // check that references are valid
  if (!is_verified())
  {
    return;
  }

  // skip if trivial
  if (gtm_bco == m_bco_reference.second)
  {
    return;
  }

  auto predicted_opt = get_predicted_fee_bco(gtm_bco);
  if (!predicted_opt)
  {
    // Safety belt: print error and exit if optional is not set
    std::cout << PHWHERE << " No predicted fee_bco for given gtm_bco, exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  const uint32_t fee_bco_predicted = predicted_opt.value();
  const double delta_fee_bco = get_signed_fee_bco_diff(fee_bco,fee_bco_predicted);
  const double gtm_bco_difference = get_signed_gtm_bco_diff(gtm_bco,m_bco_reference.second);

  m_multiplier_adjustment_numerator += gtm_bco_difference*delta_fee_bco;
  m_multiplier_adjustment_denominator += gtm_bco_difference*gtm_bco_difference;
  ++m_multiplier_adjustment_count;

  if (verbosity())
  {
    const auto default_precision{std::cout.precision()};
    std::cout << "MicromegasBcoMatchingInformation_v2::update_multiplier_adjustment -"
              << " m_multiplier_adjustment_count: " << m_multiplier_adjustment_count
              << std::setprecision(10)
              << " m_multiplier: " << get_adjusted_multiplier()
              << " adjustment: " << m_multiplier_adjustment_numerator / m_multiplier_adjustment_denominator
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

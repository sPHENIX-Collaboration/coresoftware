
/*!
 * \file MicromegasBcoMatchingInformation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief handles matching between GTM and FEE BCO clock
 */

#include "MicromegasBcoMatchingInformation.h"

#include <Event/packet.h>

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
      const bool is_hex = (o.flags()&std::ios_base::hex);
      o << "{ ";
      bool first = true;
      for (const auto& value : list)
      {
        if (!first)
        {
          o << ", ";
        }
        if( is_hex )
        { o << "0x"; }
        o << value;
        first = false;
      }
      o << " }";
    }
    return o;
  }

  // get the difference between two BCO.
  template<class T>
    inline static constexpr T get_bco_diff( const T& first, const T& second )
  { return first < second ? (second-first):(first-second); }

  //! copied from micromegas/MicromegasDefs.h, not available here
  static constexpr int m_nchannels_fee = 256;

}

// this is the clock multiplier from lvl1 to fee clock
/* todo: should replace with actual rational number for John K. */
double MicromegasBcoMatchingInformation::m_multiplier = 4.262916255;

//___________________________________________________
std::optional<uint32_t> MicromegasBcoMatchingInformation::get_predicted_fee_bco( uint64_t gtm_bco ) const
{
  // check proper initialization
  if( !m_verified ) { return std::nullopt; }

  // get gtm bco difference with proper rollover accounting
  const uint64_t gtm_bco_difference = (gtm_bco >= m_gtm_bco_first) ?
    (gtm_bco - m_gtm_bco_first):
    (gtm_bco + (1ULL<<40U) - m_gtm_bco_first);

  // convert to fee bco, and truncate to 20 bits
  const uint64_t fee_bco_predicted = m_fee_bco_first + m_multiplier*(gtm_bco_difference);
  return  uint32_t(fee_bco_predicted & 0xFFFFFU);
}

//___________________________________________________
void MicromegasBcoMatchingInformation::save_gtm_bco_information( Packet* packet )
{
  // append gtm_bco from taggers in this event to packet-specific list of available lv1_bco
  const int n_tagger = packet->lValue(0, "N_TAGGER");
  for (int t = 0; t < n_tagger; t++)
  {
    bool is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
    if (is_lvl1)
    {
      uint64_t gtm_bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));
      m_gtm_bco_list.push_back(gtm_bco);
    }
  }

  if(verbosity() && !m_gtm_bco_list.empty())
  {
    // get packet id
    const int packet_id = packet->getIdentifier();

    std::cout
      << "MicromegasBcoMatchingInformation::save_gtm_bco_information -"
      << " packet: " << packet_id
      << " n_tagger: " << n_tagger
      << std::endl;

    std::cout
      << "MicromegasBcoMatchingInformation::save_gtm_bco_information -"
      << " packet: " << packet_id
      << " gtm_bco: " << std::hex << m_gtm_bco_list << std::dec
      << std::endl;

    // also print predicted fee bco
    if( m_verified )
    {
      std::list<uint32_t> fee_bco_predicted_list;
      std::transform(
        m_gtm_bco_list.begin(),
        m_gtm_bco_list.end(),
        std::back_inserter(fee_bco_predicted_list),
        [this](const uint64_t& gtm_bco ){ return get_predicted_fee_bco(gtm_bco).value(); } );

      std::cout
        << "MicromegasBcoMatchingInformation::save_gtm_bco_information -"
        << " packet: " << packet_id
        << " fee_bco_predicted: " << std::hex << fee_bco_predicted_list << std::dec
        << std::endl;
    }
  }
}

//___________________________________________________
bool MicromegasBcoMatchingInformation::find_reference( Packet* packet )
{
  // store gtm bco and diff to previous in an array
  std::vector<uint64_t> gtm_bco_list;
  std::vector<uint64_t> gtm_bco_diff_list;
  for( const auto& gtm_bco:m_gtm_bco_list )
  {
    if( !gtm_bco_list.empty() )
    {
      // add difference to last
      // get gtm bco difference with proper rollover accounting
      const uint64_t gtm_bco_difference = (gtm_bco >= gtm_bco_list.back()) ?
        (gtm_bco -  gtm_bco_list.back()):
        (gtm_bco + (1ULL<<40U) - gtm_bco_list.back());

      gtm_bco_diff_list.push_back(m_multiplier*gtm_bco_difference);
    }

    // append to list
    gtm_bco_list.push_back( gtm_bco );
  }

  uint32_t fee_bco_prev = 0;
  bool has_fee_bco_prev = false;

  // get number of waveforms
  const int n_waveform = packet->iValue(0, "NR_WF");
  for (int iwf = 0; iwf < n_waveform; ++iwf)
  {
    const unsigned short channel = packet->iValue( iwf, "CHANNEL" );

    // bound check
    if( channel >= m_nchannels_fee )
    { continue; }

    const uint32_t fee_bco = static_cast<uint32_t>(packet->iValue(iwf, "BCO"));
    if( !has_fee_bco_prev )
    {
      fee_bco_prev = fee_bco;
      has_fee_bco_prev = true;
    }

    // calculate difference
    const uint64_t fee_bco_diff = get_bco_diff( fee_bco, fee_bco_prev );

    // discard identical fee_bco
    if( fee_bco_diff < m_max_fee_bco_diff )
    { continue; }

    // loop for matching diff in gtm_bco array
    for( size_t i = 0; i < gtm_bco_diff_list.size(); ++i )
    {
      if( get_bco_diff( gtm_bco_diff_list[i], fee_bco_diff ) < m_max_fee_bco_diff )
      {
        m_verified = true;
        m_gtm_bco_first = gtm_bco_list[i];
        m_fee_bco_first = fee_bco_prev;

        if( verbosity() )
        {
          std::cout << "MicromegasBcoMatchingInformation::find_reference - matching is verified" << std::endl;
          std::cout
            << "MicromegasBcoMatchingInformation::find_reference -"
            << " m_gtm_bco_first: " << std::hex << m_gtm_bco_first << std::dec
            << std::endl;
          std::cout
            << "MicromegasBcoMatchingInformation::find_reference -"
            << " m_fee_bco_first: " << std::hex << m_fee_bco_first << std::dec
            << std::endl;
        }
        return true;
      }
    }
    // update previous fee_bco
    fee_bco_prev = fee_bco;
  }
  return false;
}


//___________________________________________________
std::optional<uint64_t> MicromegasBcoMatchingInformation::find_gtm_bco( uint32_t fee_bco )
{
  // make sure the bco matching is properly initialized
  if( !m_verified )
  {
    return std::nullopt;
  }
  // find matching gtm bco in map
  const auto bco_matching_iter = std::find_if(
    m_bco_matching_list.begin(),
    m_bco_matching_list.end(),
    [fee_bco]( const m_bco_matching_pair_t& pair )
    { return get_bco_diff( pair.first, fee_bco ) < m_max_fee_bco_diff; } );

  if( bco_matching_iter != m_bco_matching_list.end() )
  {

    return bco_matching_iter->second;

  } else {

    // find element for which predicted fee_bco is the closest to request
    const auto iter = std::find_if(
      m_gtm_bco_list.begin(),
      m_gtm_bco_list.end(),
      [this, fee_bco]( const uint64_t& gtm_bco )
      { return get_bco_diff( get_predicted_fee_bco(gtm_bco).value(), fee_bco ) <  m_max_gtm_bco_diff; } );

    // check
    if( iter != m_gtm_bco_list.end() )
    {
      const auto gtm_bco = *iter;
      if (verbosity())
      {
        const auto predicted_fee_bco = get_predicted_fee_bco(gtm_bco).value();
        const auto fee_bco_diff = get_bco_diff(predicted_fee_bco, fee_bco);

        std::cout << "MicromegasBcoMatchingInformation::find_gtm_bco -"
          << std::hex
          << " fee_bco: 0x" << fee_bco
          << " predicted: 0x" << predicted_fee_bco
          << " gtm_bco: 0x" << gtm_bco
          << std::dec
          << " difference: " << fee_bco_diff
          << std::endl;
      }

      // save fee_bco and gtm_bco matching in map
      m_bco_matching_list.emplace_back(fee_bco, gtm_bco);

      // remove gtm bco from runing list
      m_gtm_bco_list.erase(iter);

      return gtm_bco;
    } else {

      if(verbosity() && m_orphans.insert(fee_bco).second)
      {

        // find element for which predicted fee_bco is the closest to request
        const auto iter2 = std::min_element(
          m_gtm_bco_list.begin(),
          m_gtm_bco_list.end(),
          [this, fee_bco]( const uint64_t& first, const uint64_t& second )
          { return get_bco_diff( get_predicted_fee_bco(first).value(), fee_bco ) <  get_bco_diff( get_predicted_fee_bco(second).value(), fee_bco ); } );

        const int fee_bco_diff = (iter2 != m_gtm_bco_list.end()) ?
          get_bco_diff( get_predicted_fee_bco(*iter2).value(), fee_bco ):-1;

        std::cout << "MicromegasBcoMatchingInformation::find_gtm_bco -"
          << std::hex
          << " fee_bco: 0x" << fee_bco
          << std::dec
          << " gtm_bco: none"
          << " difference: " << fee_bco_diff
          << std::endl;
      }
      return std::nullopt;
    }
  }

  // never reached
  return std::nullopt;

}

//___________________________________________________
void MicromegasBcoMatchingInformation::cleanup()
{
  // remove old gtm_bco and matching
  while( m_gtm_bco_list.size() > m_max_matching_data_size ) { m_gtm_bco_list.pop_front(); }
  while( m_bco_matching_list.size() > m_max_matching_data_size ) { m_bco_matching_list.pop_front(); }

  // clear orphans
  m_orphans.clear();
}

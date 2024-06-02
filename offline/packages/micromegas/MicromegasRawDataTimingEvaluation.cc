/*!
 * \file MicromegasRawDataTimingEvaluation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasRawDataTimingEvaluation.h"
#include "MicromegasDefs.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <cassert>
#include <fstream>
#include <memory>
#include <set>

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

  // get the difference between two BCO.
  template<class T>
    inline static constexpr T get_bco_diff( const T& first, const T& second )
  { return first < second ? (second-first):(first-second); }

}

// this is the clock multiplier from lvl1 to fee clock
/* todo: should replace with actual rational number for John K. */
double MicromegasRawDataTimingEvaluation::bco_matching_information_t::m_multiplier = 4.262916255;

//_________________________________________________________
void MicromegasRawDataTimingEvaluation::bco_matching_information_t::truncate( unsigned int maxsize )
{
  while( m_gtm_bco_list.size() > maxsize ) { m_gtm_bco_list.pop_front(); }
  while( m_bco_matching_list.size() > maxsize ) { m_bco_matching_list.pop_front(); }
}

//_________________________________________________________
unsigned int MicromegasRawDataTimingEvaluation::bco_matching_information_t::get_predicted_fee_bco( uint64_t gtm_bco ) const
{
  // check proper initialization
  if( !(m_has_gtm_bco_first && m_has_fee_bco_first ) ) { return 0; }

  // get gtm bco difference with proper rollover accounting
  uint64_t gtm_bco_difference = (gtm_bco >= m_gtm_bco_first) ?
    (gtm_bco - m_gtm_bco_first):
    (gtm_bco + (1ULL<<40U) - m_gtm_bco_first);

  // convert to fee bco, and truncate to 20 bits
  uint64_t fee_bco_predicted = m_fee_bco_first + m_multiplier*(gtm_bco_difference);
  return (unsigned int)(fee_bco_predicted & 0xFFFFFU);
}

//_________________________________________________________
void MicromegasRawDataTimingEvaluation::Container::Reset()
{ waveforms.clear(); }

//_________________________________________________________
MicromegasRawDataTimingEvaluation::MicromegasRawDataTimingEvaluation(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int MicromegasRawDataTimingEvaluation::Init(PHCompositeNode* /*topNode*/)
{

  std::cout << "MicromegasRawDataTimingEvaluation::Init -"
    << " bco_matching_information_t::m_multiplier: " << bco_matching_information_t::m_multiplier
    << std::endl;

  m_evaluation_file.reset(new TFile(m_evaluation_filename.c_str(), "RECREATE"));
  m_evaluation_tree = new TTree("T", "T");
  m_container = new Container;
  m_evaluation_tree->Branch("Event", &m_container);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasRawDataTimingEvaluation::InitRun(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int MicromegasRawDataTimingEvaluation::process_event(PHCompositeNode* topNode)
{
  // load relevant nodes
  // PRDF node
  auto event = findNode::getClass<Event>(topNode, "PRDF");
  assert(event);

  // check event type
  if (event->getEvtType() >= 8)
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  m_container->Reset();

  // loop over TPOT packets
  for (const auto& packet_id : MicromegasDefs::m_packet_ids)
  {
    std::unique_ptr<Packet> packet(event->getPacket(packet_id));
    if (!packet)
    {
      // no data
      if (Verbosity() > 1)
      {
        std::cout << "MicromegasRawDataTimingEvaluation::process_event - packet " << packet_id << " not found." << std::endl;
      }
      continue;
    }

    // get relevant bco matching information
    auto& bco_matching_information = m_bco_matching_information_map[packet_id];

    // append gtm_bco from taggers in this event to packet-specific list of available lv1_bco
    const int n_tagger = packet->lValue(0, "N_TAGGER");
    for (int t = 0; t < n_tagger; t++)
    {
      bool is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
      uint64_t gtm_bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));

      if (is_lvl1)
      {
        // initialize first gtm_bco
        if( !bco_matching_information.m_has_gtm_bco_first )
        {
          bco_matching_information.m_gtm_bco_first = gtm_bco;
          bco_matching_information.m_has_gtm_bco_first = true;
          if( Verbosity() )
          {
            std::cout
              << "MicromegasRawDataTimingEvaluation::process_event -"
              << " packet: " << packet_id
              << std::hex
              << " m_gtm_bco_first: 0x" << gtm_bco
              << std::dec
              << std::endl;
          }
        }

        // store in list of available bco
        bco_matching_information.m_gtm_bco_list.push_back(gtm_bco);
      }
    }

    // get number of waveforms
    const auto n_waveform = packet->iValue(0, "NR_WF");
    m_waveform_count_total += n_waveform;

    if (Verbosity())
    {
      std::cout << "MicromegasRawDataTimingEvaluation::process_event -"
                << " packet: " << packet_id
                << " taggers: " << n_tagger
                << " n_gtm_bco: " << bco_matching_information.m_gtm_bco_list.size()
                << " n_waveform: " << n_waveform
                << std::endl;

      if (!bco_matching_information.m_gtm_bco_list.empty())
      {
        std::cout
          << "MicromegasRawDataTimingEvaluation::process_event -"
          << " packet: " << packet_id
          << " gtm_bco: " << std::hex << bco_matching_information.m_gtm_bco_list << std::dec
          << std::endl;

        // also print predicted fee bco
        std::list<unsigned int> fee_bco_predicted_list;
        std::transform(
          bco_matching_information.m_gtm_bco_list.begin(),
          bco_matching_information.m_gtm_bco_list.end(),
          std::back_inserter(fee_bco_predicted_list),
          [&bco_matching_information](const uint64_t& gtm_bco ){ return bco_matching_information.get_predicted_fee_bco(gtm_bco); } );

        std::cout
          << "MicromegasRawDataTimingEvaluation::process_event -"
          << " packet: " << packet_id
          << " fee_bco_predicted: " << std::hex << fee_bco_predicted_list << std::dec
          << std::endl;
      }
    }

    // keep track of orphans
    using fee_pair_t = std::pair<unsigned int, unsigned int>;
    std::set<fee_pair_t> orphans;

    for (int iwf = 0; iwf < n_waveform; ++iwf)
    {
      Waveform waveform;
      waveform.packet_id = packet_id;
      waveform.fee_id = packet->iValue(iwf, "FEE");
      waveform.channel = packet->iValue( iwf, "CHANNEL" );

      // bound check
      if( waveform.channel >= MicromegasDefs::m_nchannels_fee )
      {
        if( Verbosity() )
        { std::cout << "MicromegasRawDataTimingEvaluation::process_event - invalid channel: " << waveform.channel << std::endl; }
        continue;
      }

      // beam crossing
      waveform.fee_bco = static_cast<uint32_t>(packet->iValue(iwf, "BCO"));
      waveform.gtm_bco = 0;

      // initialize first fee_bco
      if( !bco_matching_information.m_has_fee_bco_first )
      {
        bco_matching_information.m_fee_bco_first = waveform.fee_bco;
        bco_matching_information.m_has_fee_bco_first = true;
        if( Verbosity() )
        {
          std::cout
            << "MicromegasRawDataTimingEvaluation::process_event -"
            << " packet: " << packet_id
            << std::hex
            << " m_fee_bco_first: 0x" << waveform.fee_bco
            << std::dec
            << std::endl;
        }
      }

      // find matching gtm bco
      const auto bco_matching_iter = std::find_if(
        bco_matching_information.m_bco_matching_list.begin(),
        bco_matching_information.m_bco_matching_list.end(),
        [&waveform]( const m_bco_matching_pair_t& pair )
        { return get_bco_diff( pair.first, waveform.fee_bco ) < bco_matching_information_t::m_max_fee_bco_diff; } );

      if( bco_matching_iter != bco_matching_information.m_bco_matching_list.end() )
      {

        // found matching gtm
        waveform.gtm_bco = bco_matching_iter->second;
        waveform.fee_bco_predicted = bco_matching_information.get_predicted_fee_bco(bco_matching_iter->second);

      } else {

        const auto iter = std::min_element(
          bco_matching_information.m_gtm_bco_list.begin(),
          bco_matching_information.m_gtm_bco_list.end(),
          [&waveform, &bco_matching_information]( const uint64_t& first, const uint64_t& second )
          { return get_bco_diff( bco_matching_information.get_predicted_fee_bco(first), waveform.fee_bco ) <  get_bco_diff( bco_matching_information.get_predicted_fee_bco(second), waveform.fee_bco ); } );

        if( iter != bco_matching_information.m_gtm_bco_list.end() && get_bco_diff( bco_matching_information.get_predicted_fee_bco(*iter), waveform.fee_bco ) < bco_matching_information_t::m_max_gtm_bco_diff )
        {
          const auto gtm_bco = *iter;
          if (false && Verbosity())
          {
            std::cout << "MicromegasRawDataTimingEvaluation::process_event -"
              << " fee_id: " << waveform.fee_id
              << std::hex
              << " fee_bco: 0x" << waveform.fee_bco
              << " predicted: 0x" << bco_matching_information.get_predicted_fee_bco(gtm_bco)
              << " gtm_bco: 0x" << gtm_bco
              << std::dec
              << " difference: " << get_bco_diff(bco_matching_information.get_predicted_fee_bco(gtm_bco), waveform.fee_bco)
              << std::endl;
          }

          // save fee_bco and gtm_bco matching in map
          bco_matching_information.m_bco_matching_list.emplace_back(waveform.fee_bco, gtm_bco);

          // store
          waveform.gtm_bco = gtm_bco;
          waveform.fee_bco_predicted = bco_matching_information.get_predicted_fee_bco(gtm_bco);

          // remove bco from running list
          bco_matching_information.m_gtm_bco_list.erase(iter);

          /*
          * if matching information is not verified, and the found match is not trivial (0),
          * change verified flag to true.
          */
          if( !bco_matching_information.m_verified && get_bco_diff( waveform.fee_bco, bco_matching_information.m_fee_bco_first ) > bco_matching_information_t::m_max_gtm_bco_diff )
          {
            bco_matching_information.m_verified = true;
          }

        } else {
          if (Verbosity() && orphans.insert(std::make_pair(waveform.fee_id, waveform.fee_bco)).second)
          {

            const int gtm_bco_diff = (iter != bco_matching_information.m_gtm_bco_list.end()) ?
              get_bco_diff( bco_matching_information.get_predicted_fee_bco(*iter), waveform.fee_bco ):-1;

            std::cout << "MicromegasRawDataTimingEvaluation::process_event -"
              << " fee_id: " << waveform.fee_id
              << std::hex
              << " fee_bco: 0x" << waveform.fee_bco
              << std::dec
              << " gtm_bco: none"
              << " difference: " << gtm_bco_diff
              << std::endl;
          }

          // increment count
          ++m_waveform_count_dropped;

          /*
          * if no match is found, and matching_information has not been verified,
          * try using this BCO as a reference instead
          */
          if( !bco_matching_information.m_verified && waveform.fee_bco > bco_matching_information.m_fee_bco_first && waveform.fee_id != 11 )
          {
            std::cout << "MicromegasRawDataTimingEvaluation::process_event - adjusting FEE reference" << std::endl;
            bco_matching_information.m_has_fee_bco_first = true;
            bco_matching_information.m_fee_bco_first = waveform.fee_bco;
          }
        }
      }

      // increment number of waveforms found for this gtm_bco
      ++m_bco_map[waveform.gtm_bco];
      m_container->waveforms.push_back(waveform);
    }

    // cleanup
    bco_matching_information.truncate(bco_matching_information_t::m_max_matching_data_size);

  }

  // fill evaluation tree
  m_evaluation_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasRawDataTimingEvaluation::End(PHCompositeNode* /*topNode*/)
{
  if (m_evaluation_file && m_evaluation_tree)
  {
    m_evaluation_file->cd();
    m_evaluation_tree->Write();
    m_evaluation_file->Close();
  }

  // print bco map
  if (false && Verbosity())
  {
    for (const auto& [bco, nwaveforms] : m_bco_map)
    {
      std::cout
        << "MicromegasRawDataTimingEvaluation::End - "
        << " bco: 0x" << std::hex << bco << std::dec
        << ", nwaveforms: " << nwaveforms
        << std::endl;
    }
  }

  // print bco list, for offline processing
  if (false && Verbosity())
  {
    std::cout << "const std::vector<uint64_t> gtm_bco_list = {" << std::endl;
    bool first = true;
    int count = 0;
    for (const auto& [bco, nwaveforms] : m_bco_map)
    {
      if (!first)
      {
        std::cout << ", ";
      }
      first = false;
      if (count == 10)
      {
        count = 0;
        std::cout << std::endl;
      }
      std::cout << " 0x" << std::hex << bco << std::dec;
      ++count;
    }
    std::cout << std::endl
              << "};" << std::endl;
  }

  if( Verbosity() )
  {
    std::cout << "MicromegasRawDataTimingEvaluation::End - waveform_count_total: " << m_waveform_count_total << std::endl;
    std::cout << "MicromegasRawDataTimingEvaluation::End - waveform_count_dropped: " << m_waveform_count_dropped << std::endl;
    std::cout << "MicromegasRawDataTimingEvaluation::End - ratio: " << double(m_waveform_count_dropped)/m_waveform_count_total << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

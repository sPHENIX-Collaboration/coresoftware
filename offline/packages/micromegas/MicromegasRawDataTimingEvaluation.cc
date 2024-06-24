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
  {
    const auto default_precision{std::cout.precision()};
    std::cout << "MicromegasRawDataTimingEvaluation::Init -"
      << std::setprecision(10)
      << " MicromegasBcoMatchingInformation::multiplier: " << MicromegasBcoMatchingInformation::get_gtm_clock_multiplier()
      << std::setprecision(default_precision)
      << std::endl;
    }

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
    bco_matching_information.set_verbosity( Verbosity() );

    // read gtm bco information
    bco_matching_information.save_gtm_bco_information( packet.get() );

    // get number of waveforms
    const auto n_waveform = packet->iValue(0, "NR_WF");
    m_waveform_count_total[packet_id] += n_waveform;

    if (Verbosity())
    {
      std::cout
        << "MicromegasRawDataTimingEvaluation::process_event -"
        << " packet: " << packet_id
        << " n_waveform: " << n_waveform
        << std::endl;
      bco_matching_information.print_gtm_bco_information();
    }

    // try find reference
    if( !bco_matching_information.is_verified() )
    { bco_matching_information.find_reference( packet.get() ); }

    // drop packet if not found
    if( !bco_matching_information.is_verified() )
    {
      std::cout << "MicromegasRawDataTimingEvaluation::process_event - bco_matching not verified, dropping packet" << std::endl;
      m_waveform_count_dropped[packet_id] += n_waveform;
      bco_matching_information.cleanup();
      continue;
    }

    for (int iwf = 0; iwf < n_waveform; ++iwf)
    {
      Waveform waveform;
      waveform.packet_id = packet_id;
      waveform.fee_id = packet->iValue(iwf, "FEE");
      waveform.channel = packet->iValue( iwf, "CHANNEL" );
      waveform.type = packet->iValue(iwf, "TYPE");

      // bound check
      if( waveform.channel >= MicromegasDefs::m_nchannels_fee )
      {
        if( Verbosity() )
        { std::cout << "MicromegasRawDataTimingEvaluation::process_event - invalid channel: " << waveform.channel << std::endl; }
        continue;
      }

      // beam crossing
      waveform.fee_bco = static_cast<uint32_t>(packet->iValue(iwf, "BCO"));

      const auto result = bco_matching_information.find_gtm_bco( waveform.fee_bco );
      if( result )
      {
        // found matching gtm
        waveform.gtm_bco = result.value();

        // also save predicted bco
        waveform.fee_bco_predicted = bco_matching_information.get_predicted_fee_bco( waveform.gtm_bco ).value();

      } else {

        // invalid
        waveform.gtm_bco = 0;

        // increment drop count
        ++m_waveform_count_dropped[packet_id];

      }

      // increment number of waveforms found for this gtm_bco
      ++m_bco_map[waveform.gtm_bco];
      m_container->waveforms.push_back(waveform);
    }

    // cleanup
    bco_matching_information.cleanup();

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

  for( const auto& [packet,counts]:m_waveform_count_total )
  {
    const auto& dropped = m_waveform_count_dropped[packet];
    std::cout << "MicromegasRawDataTimingEvaluation::End - packet: " << packet << std::endl;
    std::cout << "MicromegasRawDataTimingEvaluation::End - waveform_count_total: " << counts << std::endl;
    std::cout << "MicromegasRawDataTimingEvaluation::End - waveform_count_dropped: " << dropped << std::endl;
    std::cout << "MicromegasRawDataTimingEvaluation::End - ratio: " << double(dropped)/counts << std::endl;
    std::cout << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

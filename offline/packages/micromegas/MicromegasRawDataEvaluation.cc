/*!
 * \file MicromegasRawDataEvaluation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasRawDataEvaluation.h"
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
void MicromegasRawDataEvaluation::Waveform::copy_from(const MicromegasRawDataEvaluation::Sample& sample)
{
  packet_id = sample.packet_id;
  gtm_bco = sample.gtm_bco;
  fee_bco = sample.fee_bco;
  fee_bco_predicted = sample.fee_bco_predicted;
  checksum = sample.checksum;
  checksum_error = sample.checksum_error;
  fee_id = sample.fee_id;
  layer = sample.layer;
  tile = sample.tile;
  sampa_address = sample.sampa_address;
  sampa_channel = sample.sampa_channel;
  channel = sample.channel;
  strip = sample.strip;
  sample_max = sample.sample;
  adc_max = sample.adc;
  pedestal = sample.pedestal;
  rms = sample.rms;
}

//_________________________________________________________
void MicromegasRawDataEvaluation::Container::Reset()
{
  samples.clear();
  waveforms.clear();
  taggers.clear();
}

//_________________________________________________________
MicromegasRawDataEvaluation::MicromegasRawDataEvaluation(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int MicromegasRawDataEvaluation::Init(PHCompositeNode* /*topNode*/)
{
  // read calibrations
  m_calibration_data.read(m_calibration_filename);

  m_evaluation_file.reset(new TFile(m_evaluation_filename.c_str(), "RECREATE"));
  m_evaluation_tree = new TTree("T", "T");
  m_container = new Container;
  m_evaluation_tree->Branch("Event", &m_container);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasRawDataEvaluation::InitRun(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int MicromegasRawDataEvaluation::process_event(PHCompositeNode* topNode)
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

  // temporary storage for samples and waveforms, sorted by gtm bco
  std::multimap<uint64_t, Sample> sample_map;
  std::multimap<uint64_t, Waveform> waveform_map;

  // loop over TPOT packets
  for (const auto& packet_id : MicromegasDefs::m_packet_ids)
  {
    std::unique_ptr<Packet> packet(event->getPacket(packet_id));
    if (!packet)
    {
      // no data
      if (Verbosity() > 1)
      {
        std::cout << "MicromegasRawDataEvaluation::process_event - packet " << packet_id << " not found." << std::endl;
      }
      continue;
    }

    // get relevant bco matching information
    auto& bco_matching_information = m_bco_matching_information_map[packet_id];
    bco_matching_information.set_verbosity( Verbosity() );

    // read gtm bco information
    bco_matching_information.save_gtm_bco_information( packet.get() );

    // append gtm_bco from taggers in this event to packet-specific list of available lv1_bco
    const int n_tagger = packet->lValue(0, "N_TAGGER");
    for (int t = 0; t < n_tagger; t++)
    {
      TaggerInformation tagger;
      tagger.packet_id = packet_id;
      tagger.tagger_type = (uint16_t) (packet->lValue(t, "TAGGER_TYPE"));
      tagger.is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
      tagger.is_endat = static_cast<uint8_t>(packet->lValue(t, "IS_ENDAT"));
      tagger.bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));
      tagger.last_bco = static_cast<uint64_t>(packet->lValue(t, "LAST_BCO"));
      tagger.lvl1_count = static_cast<uint32_t>(packet->lValue(t, "LEVEL1_COUNT"));
      tagger.endat_count = static_cast<uint32_t>(packet->lValue(t, "ENDAT_COUNT"));

      if (m_flags & EvalTagger)
      {
        m_container->taggers.push_back(tagger);
      }
    }

    // get number of waveforms
    const auto n_waveform = packet->iValue(0, "NR_WF");
    m_waveform_count_total[packet_id] += n_waveform;

    if (Verbosity())
    {
      std::cout << "MicromegasRawDataEvaluation::process_event -"
                << " packet: " << packet_id
                << " n_waveform: " << n_waveform
                << std::endl;
    }

    // find reference from modebits, using BX_COUNTER_SYNC_T
    /*
     * This needs to be done even if the bco matching information is already verified
     * because any BX_COUNTER_SYNC_T event will break past references
     */
    bco_matching_information.find_reference_from_modebits(packet.get());

    // if bco matching information is not verified, try find reference from data
    if (!bco_matching_information.is_verified())
    {
      bco_matching_information.find_reference_from_data(packet.get());
    }

    if (m_flags & (EvalSample | EvalWaveform))
    {
      // keep track of orphans
      using fee_pair_t = std::pair<unsigned int, unsigned int>;
      std::set<fee_pair_t> orphans;

      for (int iwf = 0; iwf < n_waveform; ++iwf)
      {

        const int type = packet->iValue(iwf, "TYPE");
        if( type == MicromegasDefs::HEARTBEAT_T ) continue;

        // create running sample, assign packet, fee, layer and tile id
        Sample sample;
        sample.packet_id = packet_id;

        // get fee id, apply mapping to current fiber set, for backward compatibility
        sample.fee_id = m_mapping.get_new_fee_id(packet->iValue(iwf, "FEE"));
        const auto hitsetkey = m_mapping.get_hitsetkey(sample.fee_id);
        sample.layer = TrkrDefs::getLayer( hitsetkey );
        sample.tile = MicromegasDefs::getTileId( hitsetkey );

        // get channel
        sample.channel = packet->iValue( iwf, "CHANNEL" );

        // bound check
        if( sample.channel >= MicromegasDefs::m_nchannels_fee )
        {
          if( Verbosity() )
          { std::cout << "MicromegasRawDataEvaluation::process_event - invalid channel: " << sample.channel << std::endl; }
          continue;
        }

        // beam crossing
        sample.fee_bco = static_cast<uint32_t>(packet->iValue(iwf, "BCO"));
        sample.gtm_bco = 0;

        const auto result = bco_matching_information.find_gtm_bco( sample.fee_bco );
        if(result)
        {
          sample.gtm_bco = result.value();
          sample.fee_bco_predicted = bco_matching_information.get_predicted_fee_bco(sample.gtm_bco).value();

        } else {

          // increment count
          ++m_waveform_count_dropped[packet_id];

        }

        // increment number of waveforms found for this gtm_bco
        ++m_bco_map[sample.gtm_bco];

        // channel, sampa_channel, sampa address and strip
        sample.sampa_address = packet->iValue( iwf, "SAMPAADDRESS" );
        sample.sampa_channel = packet->iValue( iwf, "SAMPACHANNEL" );
        sample.strip = m_mapping.get_physical_strip(sample.fee_id, sample.channel);

        // get channel rms and pedestal from calibration data
        const double pedestal = m_calibration_data.get_pedestal(sample.fee_id, sample.channel);
        const double rms = m_calibration_data.get_rms(sample.fee_id, sample.channel);
        sample.pedestal = pedestal;
        sample.rms = rms;

        // get number of samples and loop
        const unsigned short samples = packet->iValue(iwf, "SAMPLES");
        if (Verbosity() > 1)
        {
          std::cout << "MicromegasRawDataEvaluation::process_event -"
                    << " fee: " << sample.fee_id
                    << " tile: " << sample.tile
                    << " layer: " << sample.layer
                    << " tile: " << sample.tile
                    << " gtm_bco: " << sample.gtm_bco
                    << " fee_bco: " << sample.fee_bco
                    << " error: " << sample.checksum_error
                    << " channel: " << sample.channel
                    << " strip: " << sample.strip
                    << " samples: " << samples
                    << std::endl;
        }

        Sample sample_max;
        for (unsigned short is = 0; is < std::min<unsigned short>(samples, 1024); ++is)
        {
          // assign sample id and corresponding adc, save copy in container
          const uint16_t adc = packet->iValue(iwf, is);
          if (adc == MicromegasDefs::m_adc_invalid)
          {
            continue;
          }
          sample.sample = is;
          sample.adc = adc;
          sample_map.emplace(sample.gtm_bco, sample);

          if (sample.adc > sample_max.adc)
          {
            sample_max = sample;
          }
        }

        if (m_flags & EvalWaveform)
        {
          Waveform waveform(sample_max);

          waveform.is_signal =
              rms > 0 &&
              waveform.adc_max >= m_min_adc &&
              waveform.sample_max >= m_sample_min &&
              waveform.sample_max < m_sample_max &&
              waveform.adc_max > pedestal + m_n_sigma * rms;

          waveform_map.emplace(waveform.gtm_bco, waveform);
        }
      }
    }

    // cleanup
    bco_matching_information.cleanup();

  }

  // copy all samples and waveform to container
  if (m_flags & EvalSample)
  {
    for (auto&& [lvl_bco, sample] : sample_map)
    {
      m_container->samples.push_back(sample);
    }
  }

  if (m_flags & EvalWaveform)
  {
    for (auto&& [gtm_bco, waveform] : waveform_map)
    {
      m_container->waveforms.push_back(waveform);
    }
  }

  // fill evaluation tree
  m_evaluation_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasRawDataEvaluation::End(PHCompositeNode* /*topNode*/)
{
  if (m_evaluation_file && m_evaluation_tree)
  {
    m_evaluation_file->cd();
    m_evaluation_tree->Write();
    m_evaluation_file->Close();
  }

  // print bco map
  if (Verbosity())
  {
    for (const auto& [bco, nwaveforms] : m_bco_map)
    {
      std::cout
        << "MicromegasRawDataEvaluation::End - "
        << " bco: 0x" << std::hex << bco << std::dec
        << ", nwaveforms: " << nwaveforms
        << std::endl;
    }
  }

  // print bco list, for offline processing
  if (Verbosity())
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
    std::cout << "MicromegasRawDataEvaluation::End - packet: " << packet << std::endl;
    std::cout << "MicromegasRawDataEvaluation::End - waveform_count_total: " << counts << std::endl;
    std::cout << "MicromegasRawDataEvaluation::End - waveform_count_dropped: " << dropped << std::endl;
    std::cout << "MicromegasRawDataEvaluation::End - ratio: " << double(dropped)/counts << std::endl;
    std::cout << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

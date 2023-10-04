/*!
 * \file MicromegasCombinedDataEvaluation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCombinedDataEvaluation.h"
#include "MicromegasDefs.h"

#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <cassert>
#include <fstream>
#include <list>
#include <memory>

namespace
{

  // streamer for lists
  template< class T >
    std::ostream& operator << ( std::ostream& out, const std::list<T>& list )
  {
    if( list.empty() ) out << "{}";
    else
    {
      out << "{ ";
      bool first = true;
      for( const auto& value:list )
      {
        if( !first ) out << ", ";
        out << value;
        first = false;
      }

      out << " }";
    }

    return out;
  }

}


//_________________________________________________________
void MicromegasCombinedDataEvaluation::Waveform::copy_from( const MicromegasCombinedDataEvaluation::Sample& sample )
{
  packet_id = sample.packet_id;
  lvl1_bco = sample.lvl1_bco;
  fee_bco = sample.fee_bco;
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
void MicromegasCombinedDataEvaluation::Container::Reset()
{
  n_tagger.clear();
  n_waveform.clear();
  samples.clear();
  waveforms.clear();
  lvl1_bco_list.clear();
  lvl1_count_list.clear();
}

//_________________________________________________________
MicromegasCombinedDataEvaluation::MicromegasCombinedDataEvaluation( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasCombinedDataEvaluation::Init(PHCompositeNode* /*topNode*/ )
{
  // read calibrations
  m_calibration_data.read( m_calibration_filename );

  m_evaluation_file.reset( new TFile( m_evaluation_filename.c_str(), "RECREATE" ) );
  m_evaluation_tree = new TTree( "T", "T" );
  m_evaluation_tree->SetAutoSave( 5000 );
  m_container = new Container;
  m_evaluation_tree->Branch( "Event", &m_container );
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasCombinedDataEvaluation::InitRun(PHCompositeNode* /*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
int MicromegasCombinedDataEvaluation::process_event(PHCompositeNode *topNode)
{

  // load raw hits container  
  auto rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, m_rawhitnodename);
  assert( rawhitcontainer );

  // reset container
  m_container->Reset();

  // temporary storage for samples and waveforms, sorted by lvl1 bco
  std::multimap<uint64_t, Sample> sample_map;
  std::multimap<uint64_t, Waveform> waveform_map;

  // map number of waveforms per packet
  std::map<unsigned int, size_t> packet_waveforms;
  
  // loop over raw hits
  if( Verbosity() )
  { std::cout << "MicromegasCombinedDataEvaluation::process_event - hits: " << rawhitcontainer->get_nhits() << std::endl; }
  
  for( unsigned int ihit = 0; ihit < rawhitcontainer->get_nhits(); ++ihit )
  {
    const auto rawhit = rawhitcontainer->get_hit(ihit);
    const auto packet_id = rawhit->get_packetid();

    // make sure packet is valid
    if( std::find( std::begin(MicromegasDefs::m_packet_ids), std::end(MicromegasDefs::m_packet_ids ), packet_id) == std::end(MicromegasDefs::m_packet_ids ) )
    {
      std::cout << "MicromegasCombinedDataEvaluation::process_event - invalid packet: " << packet_id << std::endl;
      continue;
    }
    
    ++packet_waveforms[packet_id];

    // create running sample, assign packet, fee, layer and tile id
    Sample sample;
    sample.packet_id = packet_id;
    sample.fee_id = rawhit->get_fee();

    const auto hitsetkey = m_mapping.get_hitsetkey(sample.fee_id);
    sample.layer = TrkrDefs::getLayer( hitsetkey );
    sample.tile = MicromegasDefs::getTileId( hitsetkey );

    // beam crossing
    sample.fee_bco = rawhit->get_bco();
    sample.lvl1_bco = rawhit->get_gtm_bco();

    // increment bco map
    ++m_bco_map[sample.lvl1_bco];

//     // checksum and checksum error
//     sample.checksum = rawhit->get_checksum();
//     sample.checksum_error = rawhit->get_checksum_error();

    // channel, sampa_channel, sampa address and strip
    sample.sampa_address = rawhit->get_sampaaddress();
    sample.sampa_channel = rawhit->get_sampachannel();
    sample.channel = rawhit->get_channel();
    sample.strip = m_mapping.get_physical_strip(sample.fee_id, sample.channel);

    // get channel rms and pedestal from calibration data
    const double pedestal = m_calibration_data.get_pedestal( sample.fee_id, sample.channel );
    const double rms = m_calibration_data.get_rms( sample.fee_id, sample.channel );
    sample.pedestal = pedestal;
    sample.rms = rms;
    
    // get number of samples and loop
    const auto samples = rawhit->get_samples();
    if( Verbosity() > 1 )
    {
      std::cout << "MicromegasCombinedDataEvaluation::process_event -"
        << " fee: " << sample.fee_id
        << " tile: " << sample.tile
        << " layer: " << sample.layer
        << " tile: " << sample.tile
        << " lvl1_bco: " << sample.lvl1_bco
        << " fee_bco: " << sample.fee_bco
        << " error: " << sample.checksum_error
        << " channel: " << sample.channel
        << " strip: " << sample.strip
        << " samples: " << samples
        << std::endl;
    }

    Sample sample_max;
    for( unsigned short is = 0; is < std::min<unsigned short>( samples, 100 ); ++is )
    {
      // assign sample id and corresponding adc, save copy in container
      auto adc = rawhit->get_adc(is);
      sample.sample = is;
      sample.adc = adc;
      sample_map.emplace( sample.lvl1_bco, sample );
      
      if( sample.adc > sample_max.adc )
      { sample_max = sample; }
      
    }

    // create waveform
    Waveform waveform( sample_max );
    waveform.is_signal =
      rms > 0 &&
      waveform.adc_max >= m_min_adc &&
      waveform.sample_max >= m_sample_min &&
      waveform.sample_max < m_sample_max &&
      waveform.adc_max > pedestal+m_n_sigma * rms;
    
    waveform_map.emplace( waveform.lvl1_bco, waveform );
  }

  // copy all samples and waveform to container
  for( auto&& [lvl_bco, sample]:sample_map )
  { m_container->samples.push_back(std::move(sample)); }

  for( auto&& [lvl1_bco, waveform]:waveform_map )
  { m_container->waveforms.push_back(std::move(waveform)); }

  // store number of waveforms
  for( const auto& [packet_id, n_waveforms]:packet_waveforms ) 
  { m_container->n_waveform.push_back(n_waveforms); }
  
  // fill evaluation tree
  m_evaluation_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasCombinedDataEvaluation::End(PHCompositeNode* /*topNode*/ )
{
  if( m_evaluation_file && m_evaluation_tree )
  {
    m_evaluation_file->cd();
    m_evaluation_tree->Write();
    m_evaluation_file->Close();
  }

  // print bco map
  if( Verbosity() )
  for( const auto& [bco,nwaveforms]:m_bco_map )
  { std::cout << "MicromegasCombinedDataEvaluation::End - bco: 0x" << std::hex << bco << std::dec << ", nwaveforms: " << nwaveforms << std::endl; }

  // print bco list, for offline processing
  if( Verbosity() )
  {
    std::cout << "const std::vector<uint64_t> lvl1_bco_list = {" << std::endl;
    bool first = true;
    int count = 0;
    for( const auto& [bco,nwaveforms]:m_bco_map )
    {
      if( !first ) std::cout << ", ";
      first = false;
      if( count == 10 )
      {
        count = 0;
        std::cout << std::endl;
      }
      std::cout << " 0x" << std::hex << bco << std::dec;
      ++count;
    }
    std::cout << std::endl << "};" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

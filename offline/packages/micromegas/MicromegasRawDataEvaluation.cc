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

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <cassert>
#include <fstream>
#include <memory>

//_________________________________________________________
void MicromegasRawDataEvaluation::Container::Reset()
{
  n_waveforms = 0;
  samples.clear();
}

//_________________________________________________________
MicromegasRawDataEvaluation::MicromegasRawDataEvaluation( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasRawDataEvaluation::Init(PHCompositeNode* /*topNode*/ )
{ 
  m_evaluation_file.reset( new TFile( m_evaluation_filename.c_str(), "RECREATE" ) );
  m_evaluation_tree = new TTree( "T", "T" );
  m_container = new Container;
  m_evaluation_tree->Branch( "Event", &m_container );
  return Fun4AllReturnCodes::EVENT_OK; 
}

//____________________________________________________________________________..
int MicromegasRawDataEvaluation::InitRun(PHCompositeNode* /*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
int MicromegasRawDataEvaluation::process_event(PHCompositeNode *topNode)
{
  
  // map fee id to detector index in histogram
  using fee_map_t = std::map<int,int>;
  fee_map_t fee_map = {
    {5, 0},      // SEIP
    {7, 1},      // SEIZ
    {6, 2},      // SCOP
    {8, 3},      // SCOZ
    {9, 4},      // SCIP
    
    // old mapping until May 23
    /* it is ok to keep it here, to be able to process older files */
    {10, 5},     // SCIZ 
    // updated after May 23
    {23, 5},     // SCIZ
    {24, 6},     // SWIP
    {25, 7},     // SWIZ

    {11, 8},     // NEIP
    {12, 9},     // NEIZ
    {19, 10},    // NCOP
    {18, 11},    // NCOZ
    {0, 12},     // NCIP
    {1, 13},     // NCIZ
    {15, 14},    // NWIP
    {14, 15},    // NWIZ
  };                 

  // load relevant nodes
  // PRDF node
  auto event = findNode::getClass<Event>(topNode, "PRDF");
  assert( event );

  // check event type
  if(event->getEvtType() >= 8)
  { return Fun4AllReturnCodes::DISCARDEVENT; }

  m_container->Reset();
  
  // loop over TPOT packets
  for( const auto& packet_id:MicromegasDefs::m_packet_ids )
  {
    std::unique_ptr<Packet> packet( event->getPacket(packet_id) );
    if( !packet )
    {
      // no data
      std::cout << "MicromegasRawDataEvaluation::process_event - packet " << packet_id << " not found." << std::endl;
      continue;
    }
        
    // taggers
    m_container->n_tagger = packet->lValue(0, "N_TAGGER");
    
    
    // get number of datasets (also call waveforms)
    const auto n_waveforms = packet->iValue(0, "NR_WF" );
    m_container->n_waveforms = n_waveforms;
    m_container->max_fee_count = packet->iValue(0, "MAX_FEECOUNT"); 

    // if( Verbosity() )
    { 
      std::cout << "MicromegasRawDataEvaluation::process_event -"
        << " packet: " << packet_id 
        << " max_fee_count: " << m_container->max_fee_count
        << " n_tagger: " << m_container->n_tagger       
        << " n_waveforms: " << n_waveforms 
        << std::endl; 
    }
    
    for (int t = 0; t < m_container->n_tagger; t++)
    {
      const auto tagger_type = static_cast<uint16_t>(packet->lValue(t, "TAGGER_TYPE"));
      const auto is_endat = static_cast<uint8_t>(packet->lValue(t, "IS_ENDAT"));
      const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
      const auto bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));
      const auto lvl1_count = static_cast<uint32_t>(packet->lValue(t, "LEVEL1_COUNT"));
      const auto endat_count = static_cast<uint32_t>(packet->lValue(t, "ENDAT_COUNT"));
      const auto last_bco = static_cast<uint64_t>(packet->lValue(t, "LAST_BCO"));
      // const auto modebits = static_cast<uint8_t>(packet->lValue(t, "MODEBITS"));

      // only printout the is_lvl1 triggers
      // if( is_lvl1 )
      {
        std::cout << "MicromegasRawDataEvaluation::process_event -"
          << " packet: " << packet_id 
          << " tagger: " << t
          << " type: " << tagger_type
          << " is_enddat: " << (bool) (is_endat)
          << " is_lvl1: " << (bool) (is_lvl1)
          << " bco: " << bco
          << " last bco: " << last_bco
          << " endat_count: " << endat_count
          << " lvl1_count: " << lvl1_count
          << std::endl;
      }
    }
    
    for( int iwf=0; iwf<n_waveforms; ++iwf )
    {
      // create running sample, assign packet id
      Sample sample;
      sample.packet_id = packet_id;

      // beam crossing, checksum, checksum error
      sample.bco = packet->iValue(iwf, "BCO");
      sample.checksum = packet->iValue(iwf, "CHECKSUM");
      sample.checksum_error = packet->iValue(iwf, "CHECKSUMERROR");

      // increment bco map
      ++m_bco_map[sample.bco];

      // get hitsetkey, layer and tile
      sample.fee_id = packet->iValue(iwf, "FEE" );
      const auto hitsetkey = m_mapping.get_hitsetkey(sample.fee_id);
      sample.layer = TrkrDefs::getLayer( hitsetkey );
      sample.tile = MicromegasDefs::getTileId( hitsetkey );

      // channel, sampa_channel, sampa address and strip
      sample.sampa_address = packet->iValue( iwf, "SAMPAADDRESS" );
      sample.sampa_channel = packet->iValue( iwf, "SAMPACHANNEL" );
      sample.channel = packet->iValue( iwf, "CHANNEL" );
      sample.strip = m_mapping.get_physical_strip(sample.fee_id, sample.channel);
      
      // get detector index and absolute channel number
      const unsigned short det_index = fee_map[sample.fee_id];
      sample.absolute_channel = sample.channel + det_index*MicromegasDefs::m_nchannels_fee;
      
      // get number of samples and loop
      const unsigned short samples = packet->iValue( iwf, "SAMPLES" );

      if( Verbosity() )
      {
        std::cout << "MicromegasRawDataEvaluation::process_event -"
          << " bco: " << sample.bco
          << " error: " << sample.checksum_error
          << " layer: " << sample.layer
          << " tile: " << sample.tile
          << " sampa_address: " << sample.sampa_address
          << " sampa_channel: " << sample.sampa_channel
          << " channel: " << sample.channel 
          << " strip: " << sample.strip 
          << " samples: " << samples
          << std::endl;
      }
        
      for( unsigned short is = 0; is < samples; ++is )
      {
        // assign sample id and corresponding adc, save copy in container
        unsigned short adc = packet->iValue(iwf,is); 
        sample.sample = is;
        sample.adc = adc;
        m_container->samples.push_back( sample );
                
      }
    }
  }
        
  m_evaluation_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasRawDataEvaluation::End(PHCompositeNode* /*topNode*/ )
{
  if( m_evaluation_file && m_evaluation_tree )
  {
    m_evaluation_file->cd();
    m_evaluation_tree->Write();
    m_evaluation_file->Close();
  }
  
  // print bco map
  if( Verbosity() )
  {
    for( const auto& [bco,nwaveforms]:m_bco_map )
    { std::cout << "MicromegasRawDataEvaluation::End - bco: " << bco << ", nwaveforms: " << nwaveforms << std::endl; }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/*!
 * \file MicromegasRawDataCalibration.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasRawDataCalibration.h"
#include "MicromegasCalibrationData.h"
#include "MicromegasDefs.h"
#include "MicromegasMapping.h"

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

//_________________________________________________________
MicromegasRawDataCalibration::MicromegasRawDataCalibration( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasRawDataCalibration::Init(PHCompositeNode* /*topNode*/ )
{

  MicromegasMapping();
  
  // histogram evaluation
  if( m_savehistograms ) create_histograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasRawDataCalibration::InitRun(PHCompositeNode* /*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
int MicromegasRawDataCalibration::process_event(PHCompositeNode *topNode)
{
  
  // map fee id to detector index in histogram
  using fee_map_t = std::map<int,int>;
  fee_map_t fee_map = {
    {5, 0},      // SEIP
    {7, 1},      // SEIZ
    {6, 2},      // SCOP
    {8, 3},      // SCOZ
    {9, 4},      // SCIP
    {10, 5},     // SCIZ
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


  // get TPOT packet number
  /*
   * for now it is the same packet number as the TPC: 4001.
   * To be fixed at a later stage.
   * check with Martin Purschke
   */
  auto packet = event->getPacket(MicromegasDefs::m_packet_id);
  if( !packet )
  {
    // no data
    std::cout << "MicromegasRawDataCalibration::process_event - event contains no TPOT data" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // get number of datasets (also call waveforms)
  const auto n_waveforms = packet->iValue(0, "NR_WF" );
  if( Verbosity() )
  { std::cout << "MicromegasRawDataCalibration::process_event - n_waveforms: " << n_waveforms << std::endl; }
  
  for( int i=0; i<n_waveforms; ++i )
  {
    auto channel = packet->iValue( i, "CHANNEL" );
    int fee = packet->iValue(i, "FEE" );
    int samples = packet->iValue( i, "SAMPLES" );
    if( Verbosity() )
    {
      std::cout
        << "MicromegasRawDataCalibration::process_event -"
        << " waveform: " << i
        << " fee: " << fee
        << " channel: " << channel
        << " samples: " << samples
        << std::endl;
    }

    // find relevant profile histogram 
    TProfile* profile = nullptr;
    auto piter = m_profile_map.lower_bound( fee );
    if( piter == m_profile_map.end() || fee < piter->first )
    {
      // create and insert
      profile = new TProfile( Form( "h_adc_channel_%i", fee ), "ADC vs channel;channel;adc", MicromegasDefs::m_nchannels_fee, 0, MicromegasDefs::m_nchannels_fee );
      profile->SetErrorOption( "s" );
      m_profile_map.insert(  piter, std::make_pair( fee, profile ) );      
    } else profile = piter->second;

    // fill
    for( int is = std::max( m_sample_min,0 ); is < std::min( m_sample_max,samples ); ++ is )
    { profile->Fill( channel, packet->iValue(i,is) ); }
    
    // fill evaluation histograms
    if( m_savehistograms )
    {
      if( m_h_fee_id ) m_h_fee_id->Fill(fee);

      // find fee index from map
      const auto iter = fee_map.find( fee );
      if( iter == fee_map.end() )
      {
        std::cout << "MicromegasRawDataCalibration::process_event - unable to find fee " << fee << " in map" << std::endl;
      } else {
        
        const auto fee_index = iter->second;
        const auto channel_index = fee_index*MicromegasDefs::m_nchannels_fee + channel;

        // loop over samples
        if( m_h_adc_channel )
        {
          for( int is = std::max( m_sample_min,0 ); is < std::min( m_sample_max,samples ); ++ is )
          { m_h_adc_channel->Fill( channel_index, packet->iValue(i,is) ); }
        }
      }      
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasRawDataCalibration::End(PHCompositeNode* /*topNode*/ )
{
  
  // write calibration data to ouput file
  if( m_profile_map.empty() ) 
  {
    std::cout << "MicromegasRawDataCalibration::End - no data" << std::endl;
  } else {

    // create calibration data object
    MicromegasCalibrationData calibration_data;
    for( const auto& [fee, profile]:m_profile_map )
    {
      for( int i = 0; i < profile->GetNbinsX(); ++ i )
      {
        const auto pedestal = profile->GetBinContent(i+1);
        const auto rms = profile->GetBinError(i+1);
        calibration_data.set_pedestal( fee, i, pedestal );
        calibration_data.set_rms( fee, i, rms );
      }
    }
    calibration_data.write( m_calibration_filename );
  }
  
  // save evaluation histograms
  if( m_savehistograms && m_histogramfile )
  {
    // create mean and rms histograms
    auto profile = m_h_adc_channel->ProfileX("h_adc_channel_profx", 1, -1, "s" );
    auto h_pedestal = new TH1F( "h_pedestal", "pedestal vs channel;channel;pedestal (adc)", MicromegasDefs::m_nchannels_total, 0, MicromegasDefs::m_nchannels_total );
    auto h_rms = new TH1F( "h_rms", "rms vs channel;channel;RMS (adc)", MicromegasDefs::m_nchannels_total, 0, MicromegasDefs::m_nchannels_total );
    for( int i =0; i<MicromegasDefs::m_nchannels_total; ++i )
    {
      h_pedestal->SetBinContent( i+1, profile->GetBinContent(i+1) );
      h_rms->SetBinContent(i+1, profile->GetBinError(i+1) );
    }

    m_histogramfile->cd();
    m_h_fee_id->Write();
    m_h_adc_channel->Write();
    h_pedestal->Write();
    h_rms->Write();
    m_histogramfile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasRawDataCalibration::create_histograms()
{
  std::cout << "MicromegasRawDataCalibration::create_histograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  m_h_fee_id = new TH1I( "h_fee_id", "FEE id;Fee id;entries", 10, 0, 10 );
  m_h_adc_channel = new TH2I( "h_adc_channel", "ADC vs channel;channel;adc", MicromegasDefs::m_nchannels_total, 0, MicromegasDefs::m_nchannels_total, MicromegasDefs::m_max_adc, 0, MicromegasDefs::m_max_adc );
}

/*!
 * \file MicromegasRawDataCalibration.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasRawDataCalibration.h"
#include "MicromegasCalibrationData.h"
#include "MicromegasDefs.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TProfile.h>

#include <cassert>
#include <fstream>
#include <memory>

//_________________________________________________________
MicromegasRawDataCalibration::MicromegasRawDataCalibration( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasRawDataCalibration::Init(PHCompositeNode* /*topNode*/ )
{
  // histogram evaluation
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasRawDataCalibration::InitRun(PHCompositeNode* /*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
int MicromegasRawDataCalibration::process_event(PHCompositeNode *topNode)
{

  // load relevant nodes
  // PRDF node
  auto event = findNode::getClass<Event>(topNode, "PRDF");
  assert( event );

  // check event type
  if(event->getEvtType() >= 8)
  { return Fun4AllReturnCodes::DISCARDEVENT; }


  // loop over TPOT packets
  for( const auto& packet_id:MicromegasDefs::m_packet_ids )
  {
    std::unique_ptr<Packet> packet( event->getPacket(packet_id) );
    if( !packet )
    {
      // no data
      if( Verbosity() )
      { std::cout << "MicromegasRawDataCalibration::process_event - event contains no TPOT data" << std::endl; }
      continue;
    }

    // get number of datasets (also call waveforms)
    const auto n_waveforms = packet->iValue(0, "NR_WF" );
    if( Verbosity() )
    { std::cout << "MicromegasRawDataCalibration::process_event - n_waveforms: " << n_waveforms << std::endl; }

    for( int i=0; i<n_waveforms; ++i )
    {

      const int type = packet->iValue(i, "TYPE");
      if( type == MicromegasDefs::HEARTBEAT_T ) continue;

      // get fee id, apply mapping to current fiber set, for backward compatibility
      const auto channel = packet->iValue( i, "CHANNEL" );
      const int fee_id = m_mapping.get_new_fee_id(packet->iValue(i, "FEE"));
      const int samples = packet->iValue( i, "SAMPLES" );
      if( Verbosity() )
      {
        std::cout
          << "MicromegasRawDataCalibration::process_event -"
          << " waveform: " << i
          << " fee_id: " << fee_id
          << " channel: " << channel
          << " samples: " << samples
          << std::endl;
      }

      // find relevant profile histogram
      TProfile* profile = nullptr;
      auto piter = m_profile_map.lower_bound( fee_id );
      if( piter == m_profile_map.end() || fee_id < piter->first )
      {
        // create and insert
        profile = new TProfile( Form( "h_adc_channel_%i", fee_id ), "ADC vs channel;channel;adc", MicromegasDefs::m_nchannels_fee, 0, MicromegasDefs::m_nchannels_fee );
        profile->SetErrorOption( "s" );
        m_profile_map.insert(  piter, std::make_pair( fee_id, profile ) );
      } else profile = piter->second;

      // fill
      for( int is = std::max( m_sample_min,0 ); is < std::min( m_sample_max,samples ); ++ is )
      {
        const uint16_t adc =  packet->iValue(i,is);
        if( adc != MicromegasDefs::m_adc_invalid )
        { profile->Fill( channel, adc); }
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
    for( const auto& [fee_id, profile]:m_profile_map )
    {
      for( int i = 0; i < profile->GetNbinsX(); ++ i )
      {
        const auto pedestal = profile->GetBinContent(i+1);
        const auto rms = profile->GetBinError(i+1);
        calibration_data.set_pedestal( fee_id, i, pedestal );
        calibration_data.set_rms( fee_id, i, rms );
      }
    }
    calibration_data.write( m_calibration_filename );
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

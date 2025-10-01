/*!
 * \file MicromegasCombinedDataCalibration.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCombinedDataCalibration.h"
#include "MicromegasCalibrationData.h"
#include "MicromegasDefs.h"

#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TProfile.h>
#include <TH2.h>

#include <boost/format.hpp>

#include <cassert>
#include <fstream>
#include <memory>

//_________________________________________________________
MicromegasCombinedDataCalibration::MicromegasCombinedDataCalibration( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasCombinedDataCalibration::Init(PHCompositeNode* /*topNode*/ )
{
  // histogram evaluation
  if( m_do_evaluation )
  {
    m_evaluation_file.reset( TFile::Open( m_evaluation_filename.c_str() ) );
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasCombinedDataCalibration::InitRun(PHCompositeNode* /*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//___________________________________________________________________________
int MicromegasCombinedDataCalibration::process_event(PHCompositeNode *topNode)
{

  // load raw hits container
  auto rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, m_rawhitnodename);
  assert(rawhitcontainer);

  // loop over TPOT packets
  for (unsigned int ihit = 0; ihit < rawhitcontainer->get_nhits(); ++ihit)
  {
    const auto rawhit = rawhitcontainer->get_hit(ihit);
    const auto fee_id = rawhit->get_fee();
    const auto channel = rawhit->get_channel();
    const auto sample_range = std::make_pair(rawhit->get_sample_begin(), rawhit->get_sample_end());
    if( Verbosity() )
    {
      std::cout
        << "MicromegasCombinedDataCalibration::process_event -"
        << " waveform: " << ihit
        << " fee_id: " << fee_id
        << " channel: " << channel
        << " samples: (" << sample_range.first << "," << sample_range.second << ")"
        << std::endl;
    }

    // find relevant profile histogram
    TProfile* profile = nullptr;
    auto piter = m_profile_map.lower_bound( fee_id );
    if( piter == m_profile_map.end() || fee_id < piter->first )
    {
      // create and insert
      const auto pname = (boost::format("p_adc_channel_%i") % fee_id ).str();
      profile = new TProfile( pname.c_str(), "ADC vs channel;channel;adc", MicromegasDefs::m_nchannels_fee, 0, MicromegasDefs::m_nchannels_fee );
      profile->SetErrorOption( "s" );
      m_profile_map.insert(  piter, std::make_pair( fee_id, profile ) );
    } else {
      profile = piter->second;
    }

    // find relevant 2D histogram
    TH2* histogram = nullptr;
    if( m_do_evaluation )
    {
      auto hiter = m_histogram_map.lower_bound( fee_id );
      if( hiter == m_histogram_map.end() || fee_id < hiter->first )
      {
        static constexpr int max_adc = 1100;

        // create and insert
        const auto hname = (boost::format("h_adc_channel_%i") % fee_id ).str();
        histogram = new TH2F( hname.c_str(), "ADC vs channel;channel;adc",
          MicromegasDefs::m_nchannels_fee, 0, MicromegasDefs::m_nchannels_fee,
          max_adc, 0, max_adc );
        m_histogram_map.insert(  hiter, std::make_pair( fee_id, histogram ) );
      } else {
        histogram = hiter->second;
      }
    }

    // fill
    for( auto is = std::max(m_sample_min,sample_range.first); is < std::min(m_sample_max,sample_range.second); ++ is )
    {
        const uint16_t adc =  rawhit->get_adc(is);
        if( adc != MicromegasDefs::m_adc_invalid )
        {
          if( profile ) { profile->Fill( channel, adc); }
          if( histogram ) { histogram->Fill( channel, adc); }
        }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasCombinedDataCalibration::End(PHCompositeNode* /*topNode*/ )
{

  // write calibration data to ouput file
  if( m_profile_map.empty() )
  {
    std::cout << "MicromegasCombinedDataCalibration::End - no data" << std::endl;
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

  if( m_do_evaluation && m_evaluation_file )
  {
    m_evaluation_file->cd();
    for( const auto& [feeid,h]:m_histogram_map ) { h->Write(); }
    for( const auto& [feeid,p]:m_profile_map ) { p->Write(); }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

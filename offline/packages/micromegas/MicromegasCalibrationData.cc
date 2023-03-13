/*!
 * \file MicromegasCalibrationData.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"
#include "MicromegasMapping.h"

#include <cdbobjects/CDBTTree.h>

#include <TFile.h>
#include <fstream>
#include <sstream>
#include <iostream>

namespace
{
  static const std::string m_pedestal_key = "pedestal";
  static const std::string m_rms_key = "rms";
}

//________________________________________________________________________-
void MicromegasCalibrationData::read( const std::string& filename )
{
  std::cout << "MicromegasCalibrationData::read - filename: " << filename << std::endl;

  // clear existing data
  m_calibration_map.clear();
  
  // use generic CDBTree to load 
  CDBTTree cdbttree( filename );
  cdbttree.LoadCalibrations();
    
  // loop over registered fee ids
  MicromegasMapping mapping;
  for( const auto& fee:mapping.get_fee_id_list() )
  {    
    // loop over channels
    for( int i = 0; i < m_nchannels_fee; ++i )
    {
      // read pedestal and rms
      int channel = fee*m_nchannels_fee + i;
      double pedestal = cdbttree.GetDoubleValue( channel, m_pedestal_key, 0 );
      double rms = cdbttree.GetDoubleValue( channel, m_rms_key, 0 );

      // insert in local structure
      if( !std::isnan( rms ) )
      { m_calibration_map[fee].at(i) = calibration_data_t( pedestal, rms ); }
    }
  }
}

//________________________________________________________________________-
void MicromegasCalibrationData::set_pedestal( int fee, int channel, double value )
{ m_calibration_map[fee].at(channel).m_pedestal = value; }

//________________________________________________________________________-
void MicromegasCalibrationData::set_rms( int fee, int channel, double value )
{ m_calibration_map[fee].at(channel).m_rms = value; }

//________________________________________________________________________-
void MicromegasCalibrationData::write( const std::string& filename ) const
{
  std::cout << "MicromegasCalibrationData::write - filename: " << filename << std::endl;
  if( m_calibration_map.empty() ) return;

  // use generic CDBTree to load 
  CDBTTree cdbttree( filename );
  for( const auto& [fee,array]:m_calibration_map )
  {
    for( size_t i = 0; i < array.size(); ++i )
    { 
      int channel = fee*m_nchannels_fee + i;
      const auto& pedestal =  array[i].m_pedestal;
      const auto& rms =  array[i].m_rms;
      cdbttree.SetDoubleValue( channel, m_pedestal_key, pedestal );
      cdbttree.SetDoubleValue( channel, m_rms_key, rms );
    }
  }
  
  // commit and write
  cdbttree.Commit();
  cdbttree.WriteCDBTTree();  
}

//________________________________________________________________________-
double MicromegasCalibrationData::get_pedestal( int fee, int channel ) const
{
  const auto iter = m_calibration_map.find(fee);
  return (iter != m_calibration_map.end()) ? iter->second.at(channel).m_pedestal: -1;
}

//________________________________________________________________________-
double MicromegasCalibrationData::get_rms( int fee, int channel ) const
{
  const auto iter = m_calibration_map.find(fee);
  return (iter != m_calibration_map.end()) ? iter->second.at(channel).m_rms: -1;
}

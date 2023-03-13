/*!
 * \file MicromegasCalibrationData.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"

#include <fstream>
#include <sstream>
#include <iostream>

//________________________________________________________________________-
void MicromegasCalibrationData::read( const std::string& filename )
{
  // clear existing data
  m_calibration_map.clear();

  // read file
  std::cout << "MicromegasCalibrationData::read - filename: " << filename << std::endl;
  std::ifstream in( filename.c_str() );
  for( std::string line; std::getline( in, line ); )
  {
    
    // skip commented lines
    if( !line.rfind("//", 0) ) continue;
      
    // parse line
    std::istringstream line_str( line.c_str() );
    
    int fee = 0;
    int channel = 0;
    calibration_data_t calib_data;
    line_str >> fee >> channel >> calib_data.m_pedestal >> calib_data.m_rms;
    
    // validity check
    if( !line_str || calib_data.m_rms <= 0 ) continue;
    
    // store
    m_calibration_map[fee].at(channel) = std::move(calib_data);
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
  if( m_calibration_map.empty() ) return;
  std::ofstream out( filename.c_str() );
  out << "// fee, channel, pedestal, rms" << std::endl;
  for( const auto& [fee,array]:m_calibration_map )
  {
    for( size_t i = 0; i < array.size(); ++i )
    { out << fee << " " << i << " " << array[i].m_pedestal << " " << array[i].m_rms << std::endl; }
  }
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

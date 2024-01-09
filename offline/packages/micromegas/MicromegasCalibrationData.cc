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
std::ostream& operator << (std::ostream& out, const MicromegasCalibrationData& calib_data )
{
  out << "MicromegasCalibrationData" << std::endl;
  size_t total_entries = 0;
  for( const auto& [fee_id, data_array] : calib_data.m_raw_calibration_map )
  {
    total_entries += data_array.size();
    out << "fee_id: " << fee_id << " entries: " << data_array.size() << std::endl;
    for( size_t i=0; i<data_array.size(); ++i )
    {
      const auto& data = data_array[i];
      out << "fee_id: " << fee_id << " channel: " << i << " pedestal: " << data.m_pedestal << " rms: " << data.m_rms << std::endl;
    }

  }
  out << "total entries: " << total_entries << std::endl;


  return out;
}

//________________________________________________________________________-
void MicromegasCalibrationData::read( const std::string& filename )
{
  std::cout << "MicromegasCalibrationData::read - filename: " << filename << std::endl;

  // clear existing data
  m_raw_calibration_map.clear();
  m_mapped_calibration_map.clear();

  // make sure file exists before loading, otherwise crashes
  if( !std::ifstream( filename.c_str() ).good() )
  {
    std::cout << "MicromegasCalibrationData::read -"
      << " filename: " << filename << " does not exist."
      << " No calibration loaded" << std::endl;
    return;
  }

  // use generic CDBTree to load
  CDBTTree cdbttree( filename );
  cdbttree.LoadCalibrations();

  // loop over registered fee ids
  MicromegasMapping mapping;
  for( const auto& fee:mapping.get_fee_id_list() )
  {
    // get corresponding hitset key
    const auto hitsetkey = mapping.get_hitsetkey(fee);

    // loop over channels
    for( int i = 0; i < m_nchannels_fee; ++i )
    {
      // get matching strip id
      const auto strip = mapping.get_physical_strip( fee, i );

      // read pedestal and rms
      int channel = fee*m_nchannels_fee + i;
      double pedestal = cdbttree.GetDoubleValue( channel, m_pedestal_key, 0 );
      double rms = cdbttree.GetDoubleValue( channel, m_rms_key, 0 );

      // insert in local structure
      if( !std::isnan( rms ) )
      {
        // create calibration data
        const calibration_data_t calibration_data( pedestal, rms );

        // insert in fee,channel structure
        m_raw_calibration_map[fee].at(i) = calibration_data;

        // also insert in hitsetkey, strip structure
        if( hitsetkey > 0 && strip >= 0 )
        { m_mapped_calibration_map[hitsetkey].at(strip) = calibration_data; }
      }
    }
  }
}

//________________________________________________________________________-
void MicromegasCalibrationData::set_pedestal( int fee, int channel, double value )
{ m_raw_calibration_map[fee].at(channel).m_pedestal = value; }

//________________________________________________________________________-
void MicromegasCalibrationData::set_rms( int fee, int channel, double value )
{ m_raw_calibration_map[fee].at(channel).m_rms = value; }

//________________________________________________________________________-
void MicromegasCalibrationData::write( const std::string& filename ) const
{
  std::cout << "MicromegasCalibrationData::write - filename: " << filename << std::endl;
  if( m_raw_calibration_map.empty() ) return;

  // use generic CDBTree to load
  CDBTTree cdbttree( filename );
  for( const auto& [fee,array]:m_raw_calibration_map )
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
  const auto iter = m_raw_calibration_map.find(fee);
  return (iter != m_raw_calibration_map.end()) ? iter->second.at(channel).m_pedestal: -1;
}

//________________________________________________________________________-
double MicromegasCalibrationData::get_rms( int fee, int channel ) const
{
  const auto iter = m_raw_calibration_map.find(fee);
  return (iter != m_raw_calibration_map.end()) ? iter->second.at(channel).m_rms: -1;
}

//________________________________________________________________________-
double MicromegasCalibrationData::get_pedestal_mapped( TrkrDefs::hitsetkey hitsetkey, int strip ) const
{
  const auto iter = m_mapped_calibration_map.find(hitsetkey);
  return (iter != m_mapped_calibration_map.end()) ? iter->second.at(strip).m_pedestal: -1;
}

//________________________________________________________________________-
double MicromegasCalibrationData::get_rms_mapped( TrkrDefs::hitsetkey hitsetkey, int strip ) const
{
  const auto iter = m_mapped_calibration_map.find(hitsetkey);
  return (iter != m_mapped_calibration_map.end()) ? iter->second.at(strip).m_rms: -1;
}

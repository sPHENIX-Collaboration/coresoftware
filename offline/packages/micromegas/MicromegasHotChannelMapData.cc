/*!
 * \file MicromegasHotChannelMapData.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasHotChannelMapData.h"

#include <cdbobjects/CDBTTree.h>

#include <TFile.h>
#include <fstream>
#include <sstream>
#include <iostream>

namespace
{
  static const std::string m_total_entries_key = "total_entries";
  static const std::string m_fee_id_key = "fee_id";
  static const std::string m_channel_id_key = "channel_id";
  static constexpr int m_nchannels_fee = 256;
}

//________________________________________________________________________-
std::ostream& operator << (std::ostream& out, const MicromegasHotChannelMapData& calib_data )
{
  out << "MicromegasHotChannelMapData" << std::endl;

  for( const auto& channel_id:calib_data.m_hot_channel_map )
  { out << "  " << channel_id << std::endl; }

  out << "total entries: " << calib_data.m_hot_channel_map.size() << std::endl;
  return out;
}

//________________________________________________________________________-
void MicromegasHotChannelMapData::read( const std::string& filename )
{
  std::cout << "MicromegasHotChannelMapData::read - filename: " << filename << std::endl;

  // clear existing data
  m_hot_channel_map.clear();

  // make sure file exists before loading, otherwise crashes
  if( !std::ifstream( filename.c_str() ).good() )
  {
    std::cout << "MicromegasHotChannelMapData::read -"
      << " filename: " << filename << " does not exist."
      << " No calibration loaded" << std::endl;
    return;
  }

  // use generic CDBTree to load
  CDBTTree cdbttree( filename );
  cdbttree.LoadCalibrations();

  // read total number of hot channels
  const int m_total_entries = cdbttree.GetSingleIntValue( m_total_entries_key );
  for( int i = 0; i < m_total_entries; ++ i )
  {

    // read channel id
    const int fee_id = cdbttree.GetIntValue( i, m_fee_id_key );
    const int channel_id = cdbttree.GetIntValue( i, m_channel_id_key );
    if( std::isnan(fee_id) || std::isnan( channel_id ) ) continue;

    m_hot_channel_map.emplace( fee_id, channel_id );
  }
}

//________________________________________________________________________-
void MicromegasHotChannelMapData::add_hot_channel( int fee, int channel )
{ m_hot_channel_map.emplace( fee, channel ); }

//________________________________________________________________________-
void MicromegasHotChannelMapData::write( const std::string& filename ) const
{
  std::cout << "MicromegasHotChannelMapData::write - filename: " << filename << std::endl;
  if( m_hot_channel_map.empty() ) return;

  // use generic CDBTree to load
  CDBTTree cdbttree( filename );
  cdbttree.SetSingleIntValue( m_total_entries_key, m_hot_channel_map.size() );

  int index = 0;
  for( const auto& channel_id:m_hot_channel_map )
  {
    cdbttree.SetIntValue( index, m_fee_id_key, channel_id.m_fee );
    cdbttree.SetIntValue( index, m_channel_id_key, channel_id.m_channel );
    ++index;
  }

  // commit and write
  cdbttree.Commit();
  cdbttree.WriteCDBTTree();
}

//________________________________________________________________________-
bool MicromegasHotChannelMapData::is_hot_channel( int fee, int channel ) const
{ return m_hot_channel_map.find( {fee, channel } ) != m_hot_channel_map.end(); }

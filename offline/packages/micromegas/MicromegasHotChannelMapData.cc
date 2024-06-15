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
  static const std::string m_layer_id_key = "layer_id";
  static const std::string m_tile_id_key = "tile_id";
  static const std::string m_strip_id_key = "strip_id";
}

//________________________________________________________________________-
std::ostream& operator << (std::ostream& out, const MicromegasHotChannelMapData& calib_data )
{
  out << "MicromegasHotChannelMapData" << std::endl;

  out << "total entries: " << calib_data.m_hot_channel_map.size() << std::endl;
  for( const auto& channel_id:calib_data.m_hot_channel_map )
  { out << "  " << channel_id << std::endl; }

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
    const int layer_id = cdbttree.GetIntValue( i, m_layer_id_key );
    const int tile_id = cdbttree.GetIntValue( i, m_tile_id_key );
    const int strip_id = cdbttree.GetIntValue( i, m_strip_id_key );
    if( std::isnan(layer_id) || std::isnan(tile_id) || std::isnan(strip_id) )
    { continue; }

    m_hot_channel_map.emplace( layer_id, tile_id, strip_id );
  }

  std::cout << "MicromegasHotChannelMapData::read - total entries: " << m_hot_channel_map.size() << std::endl;

}

//________________________________________________________________________-
void MicromegasHotChannelMapData::add_hot_channel( int layer, int tile, int strip )
{ m_hot_channel_map.emplace( layer, tile, strip ); }

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
    cdbttree.SetIntValue( index, m_layer_id_key, channel_id.m_layer );
    cdbttree.SetIntValue( index, m_tile_id_key, channel_id.m_tile );
    cdbttree.SetIntValue( index, m_strip_id_key, channel_id.m_strip );
    ++index;
  }

  // commit and write
  cdbttree.Commit();
  cdbttree.CommitSingle();
  cdbttree.WriteCDBTTree();
}

//________________________________________________________________________-
bool MicromegasHotChannelMapData::is_hot_channel( int layer, int tile, int strip ) const
{ return m_hot_channel_map.find({layer, tile, strip}) != m_hot_channel_map.end(); }

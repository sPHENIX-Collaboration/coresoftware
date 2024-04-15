#ifndef MICROMEGAS_MICROMEGASHOTCHANNELMAPDATA_H
#define MICROMEGAS_MICROMEGASHOTCHANNELMAPDATA_H

/*!
 * \file MicromegasHotChannelMapData.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <trackbase/TrkrDefs.h>

#include <array>
#include <map>
#include <set>
#include <string>

/// micromegas calibration data
class MicromegasHotChannelMapData
{
  public:

  /// constructor
  MicromegasHotChannelMapData() = default;

  ///@name modifiers
  //@{

  /// read calibration from file
  void read( const std::string& /*filename*/ );

  /// clear hot channels
  void clear()
  { m_hot_channel_map.clear(); }

  /// add hot channel
  void add_hot_channel( int /*layer*/, int /*tile*/, int /*strip*/ );

  //@}

  //!@name accessors
  //@{

  /// write calibration to file
  void write( const std::string& /*filename*/ ) const;

  /// returns true if channel is hot
  bool is_hot_channel( int /*layer*/, int /*tile*/, int /*strip*/ ) const;

  //@}

  private:

  /// channel id
  class channel_id_t
  {
    public:

    /// constructor
    channel_id_t( int layer, int tile, int strip ):
      m_layer(layer),
      m_tile(tile),
      m_strip(strip)
    {}

    int m_layer = 0;
    int m_tile = 0;
    int m_strip = 0;

    /// less than operator
    bool operator < (const channel_id_t& other ) const
    {
      if( m_layer != other.m_layer ) return m_layer < other.m_layer;
      else if( m_tile != other.m_tile ) return m_tile < other.m_tile;
      else return m_strip < other.m_strip;
    }

    /// streamer
    friend std::ostream& operator << ( std::ostream& out, const channel_id_t& channel_id )
    {
      out << "{ " << channel_id.m_layer << ", " << channel_id.m_tile << ", " << channel_id.m_strip << " }";
      return out;
    }

  };

  using channel_id_set_t = std::set<channel_id_t>;
  channel_id_set_t m_hot_channel_map;

  friend std::ostream& operator << (std::ostream&, const MicromegasHotChannelMapData& );

};

#endif

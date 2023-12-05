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

  /// add hot channel
  void add_hot_channel( int /*fee*/, int /*channel*/ );

  //@}

  //!@name accessors
  //@{

  /// write calibration to file
  void write( const std::string& /*filename*/ ) const;

  /// returns true if channel is hot
  bool is_hot_channel( int /*fee*/, int /*channel*/ ) const;

  //@}

  private:

  /// channel id
  class channel_id_t
  {
    public:

    /// constructor
    channel_id_t( int fee, int channel ):
      m_fee(fee),
      m_channel(channel)
    {}

    unsigned int m_fee = 0;
    unsigned int m_channel = 0;

    /// less than operator
    bool operator < (const channel_id_t& other ) const
    {
      if( m_fee != other.m_fee ) return m_fee < other.m_fee;
      else return m_channel < other.m_channel;
    }

    // streamer
    friend std::ostream& operator << ( std::ostream& out, const channel_id_t& channel_id )
    {
      out << "{ " << channel_id.m_fee << ", " << channel_id.m_channel << " }";
      return out;
    }

  };

  using channel_id_set_t = std::set<channel_id_t>;
  channel_id_set_t m_hot_channel_map;

  friend std::ostream& operator << (std::ostream&, const MicromegasHotChannelMapData& );

};

#endif

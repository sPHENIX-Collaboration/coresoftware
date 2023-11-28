#ifndef MICROMEGAS_MICROMEGASCOMBINEDDATADECODER_H
#define MICROMEGAS_MICROMEGASCOMBINEDDATADECODER_H

/*!
 * \file MicromegasCombinedDataDecoder.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"
#include "MicromegasMapping.h"

#include <fun4all/SubsysReco.h>

#include <memory>
#include <set>
#include <string>

class PHCompositeNode;

/// micromegas raw data decoder
class MicromegasCombinedDataDecoder : public SubsysReco
{
  public:

  /// constructor
  MicromegasCombinedDataDecoder( const std::string &name = "MicromegasCombinedDataDecoder" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// calibration file
  void set_calibration_file( const std::string& value ) { m_calibration_filename = value; }

  /// set number of RMS sigma used to defined static threshold on a given channel
  void set_n_sigma( double value ) { m_n_sigma = value; }

  /// set minimum ADC value, disregarding pedestal and RMS.
  /** This removes faulty channels for which calibration has failed */
  void set_min_adc( double value ) { m_min_adc = value; }

  /// set min sample for noise estimation
  void set_sample_min( int value ) { m_sample_min = value; }

  /// set min sample for noise estimation
  void set_sample_max( int value ) { m_sample_max = value; }

  /// add a hot channel to the list
  /* hot channels are defined by layer, tile and physical strip */
  void add_hot_channel( int layer, int tile, int strip )
  { m_hot_channels.emplace( layer, tile, strip ); }

  private:

  //! raw node
  std::string m_rawhitnodename = "MICROMEGASRAWHIT";

  //! calibration filename
  std::string m_calibration_filename = "TPOT_Pedestal_000.root";

  //! calibration data
  MicromegasCalibrationData m_calibration_data;

  //! mapping
  MicromegasMapping m_mapping;

  /// number of RMS sigma used to define threshold
  double m_n_sigma = 5;

  //! minimum ADC value, disregarding pedestal and RMS.
  /* This removes faulty channels for which calibration has failed */
  double m_min_adc = 50;

  /// min sample for signal
  int m_sample_min = 0;

  /// max sample for signal
  int m_sample_max = 100;

  class channel_id_t
  {
    public:

    channel_id_t( int layer, int tile, int channel ):
      _layer( layer ),
      _tile( tile ),
      _channel( channel )
    {}

    int _layer = 0;
    int _tile = 0;
    int _channel = 0;

    bool operator < (const channel_id_t& other ) const
    {
      if( _layer != other._layer ) return _layer < other._layer;
      if( _tile != other._tile ) return _tile < other._tile;
      if( _channel != other._channel ) return _channel < other._channel;
      return false;
    }

    using set = std::set<channel_id_t>;
  };

  /// hot channels
  channel_id_t::set m_hot_channels;

  /// keep track of number of hits per hitsetid
  using hitcountmap_t = std::map<TrkrDefs::hitsetkey,int>;
  hitcountmap_t m_hitcounts;

};

#endif

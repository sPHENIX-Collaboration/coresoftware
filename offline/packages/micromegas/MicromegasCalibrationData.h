#ifndef MICROMEGAS_MICROMEGASCALIBRATIONDATA_H
#define MICROMEGAS_MICROMEGASCALIBRATIONDATA_H

/*!
 * \file MicromegasCalibrationData.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <trackbase/TrkrDefs.h>

#include <array>
#include <map>
#include <string>

/// micromegas calibration data
class MicromegasCalibrationData
{
  public:

  /// constructor
  MicromegasCalibrationData() = default;

  ///@name modifiers
  //@{
  
  /// read calibration from file
  void read( const std::string& /*filename*/ );

  /// set pedestal for a given channel
  void set_pedestal( int /*fee*/, int /*channel*/, double /*value*/ );

  /// set rms for a given channel
  void set_rms( int /*fee*/, int /*channel*/, double /*value*/ );
  
  //@}
  
  //!@name accessors
  //@{
  
  /// write calibration to file
  void write( const std::string& /*filename*/ ) const;
    
  /// get pedestal for a given channel
  double get_pedestal( int /*fee*/, int /*channel*/ ) const;
 
  /// get rms for a given channel
  double get_rms( int /*fee*/, int /*channel*/ ) const;
  
  //@}
   
  private:
   
  /// simple structure to store calibration data
  class calibration_data_t
  {
    public:    
    
    calibration_data_t() = default;
    
    calibration_data_t( double pedestal, double rms ):
      m_pedestal( pedestal ),
      m_rms( rms )
    {}
    
    double m_pedestal = 0;
    double m_rms = 0;
  };

  static constexpr int m_nchannels_fee = 256;
  using calibration_vector_t = std::array<calibration_data_t,m_nchannels_fee>;

  /// map fee id to calibration vector
  using calibration_map_t = std::map<int, calibration_vector_t>;
  calibration_map_t m_calibration_map;
    
  friend std::ostream& operator << (std::ostream&, const MicromegasCalibrationData& );
  
};

#endif

#ifndef MICROMEGAS_MICROMEGASRAWDATADECODER_H
#define MICROMEGAS_MICROMEGASRAWDATADECODER_H

/*!
 * \file MicromegasRawDataDecoder.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"
#include "MicromegasMapping.h"

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;

/// micromegas raw data decoder
class MicromegasRawDataDecoder : public SubsysReco
{
  public:

  /// constructor
  MicromegasRawDataDecoder( const std::string &name = "MicromegasRawDataDecoder" );

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
  
  /// set min sample for noise estimation
  void set_sample_min( int value ) { m_sample_min = value; }
  
  /// set min sample for noise estimation
  void set_sample_max( int value ) { m_sample_max = value; }
  
  private:
  
  //! calibration filename
  std::string m_calibration_filename = "TPOT_Pedestal_000.root";
    
  //! calibration data
  MicromegasCalibrationData m_calibration_data;
 
  //! mapping
  MicromegasMapping m_mapping;
  
  /// number of RMS sigma used to define threshold
  double m_n_sigma = 5;
  
  /// min sample for signal 
  int m_sample_min = 0;
  
  /// max sample for signal
  int m_sample_max = 100;
  
  /// keep track of number of hits per hitsetid
  using hitcountmap_t = std::map<TrkrDefs::hitsetkey,int>;
  hitcountmap_t m_hitcounts;
  
};

#endif

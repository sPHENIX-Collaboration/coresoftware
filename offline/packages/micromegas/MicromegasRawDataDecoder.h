#ifndef MICROMEGAS_MICROMEGASRAWDATADECODER_H
#define MICROMEGAS_MICROMEGASRAWDATADECODER_H

/*!
 * \file MicromegasRawDataDecoder.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"

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
  
  /// calibration file
  void set_calibration_file( const std::string& value ) { m_calibration_filename = value; }
  
  private:
  
  //! calibration filename
  std::string m_calibration_filename = "TPOT_Pedestal_000.txt";
    
  //! calibration data
  MicromegasCalibrationData m_calibration_data;
  
};

#endif

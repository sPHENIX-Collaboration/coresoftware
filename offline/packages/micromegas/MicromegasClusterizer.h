#ifndef MICROMEGAS_MICROMEGASCLUSTERIZER_H
#define MICROMEGAS_MICROMEGASCLUSTERIZER_H

/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

//! micromegas clusterizer
class MicromegasClusterizer : public SubsysReco
{
 public:

  //! constructor
  MicromegasClusterizer( const std::string &name = "MicromegasClusterizer" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  /// set default pedestal
  void set_default_pedestal( double value )
  { m_default_pedestal = value; }

  /// set whether default pedestal is used or not
  void set_use_default_pedestal( bool value )
  { m_use_default_pedestal = value; }

  /// calibration file
  void set_calibration_file( const std::string& value )
  { m_calibration_filename = value; }

  private:

  //!@name calibration filename
  //@{

  /// if true, use default pedestal to get hit charge. Relies on calibration data otherwise
  bool m_use_default_pedestal = true;

  /// default pedestal
  double m_default_pedestal = 74.6;

  /// calibration filename
  std::string m_calibration_filename;

  /// calibration data
  MicromegasCalibrationData m_calibration_data;

  //@}


};

#endif

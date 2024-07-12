#ifndef MICROMEGAS_MICROMEGASCLUSTERIZER_H
#define MICROMEGAS_MICROMEGASCLUSTERIZER_H

/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>
#include <memory>
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

  void set_drop_single_strips(bool drop)
  { m_drop_single_strips = drop; }

  /// calibration file
  void set_calibration_file( const std::string& value )
  { m_calibration_filename = value; }

  /// smearing (rphi)
  void set_added_smear_sigma_rphi( double value )
  { m_added_smear_sigma_rphi = value; }

  /// smearing (z)
  void set_added_smear_sigma_z( double value )
  { m_added_smear_sigma_z = value; }

  private:

  //!@name calibration filename
  //@{

  // discard single strip clusters if true
  bool m_drop_single_strips = false;

  /// if true, use default pedestal to get hit charge. Relies on calibration data otherwise
  bool m_use_default_pedestal = true;

  /// default pedestal
  double m_default_pedestal = 74.6;

  /// calibration filename
  std::string m_calibration_filename;

  /// calibration data
  MicromegasCalibrationData m_calibration_data;

  //@}

  //! additional smearing of primary electrons (cm)
  /** it is used to adjust the Micromegas resolution to actual measurements */
  double m_added_smear_sigma_z = 0;
  double m_added_smear_sigma_rphi = 0;

  //! rng de-allocator
  class Deleter
  {
   public:
    //! deletion operator
    void operator()(gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  //! random generator that conform with sPHENIX standard
  /*! using a unique_ptr with custom Deleter ensures that the structure is properly freed when parent object is destroyed */
  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif

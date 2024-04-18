#ifndef TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
#define TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
/**
 * \file TpcSpaceChargeMatrixInversion.h
 * \brief aggregate space charge reconstruction matrices from several jobs and perform the inversion
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/Fun4AllBase.h>

#include <memory>

// forward declaration
class TpcSpaceChargeMatrixContainer;
class TpcDistortionCorrectionContainer;

/**
 * \class TpcSpaceChargeMatrixInversion
 * \brief performs space charge distortion reconstruction using tracks
 */

class TpcSpaceChargeMatrixInversion: public Fun4AllBase
{
  public:

  /// constructor
  TpcSpaceChargeMatrixInversion( const std::string& = "TPCSPACECHARGEMATRIXINVERSION" );

  ///@name modifiers
  //@{

  /// set central membrane distortion file
  void set_central_membrane_correction_filename( const std::string& );

  /// add space charge correction matrix to current. Returns true on success
  bool add( const TpcSpaceChargeMatrixContainer& );

  /// add space charge correction matrix, loaded from file, to current. Returns true on success
  bool add_from_file( const std::string& /*filename*/, const std::string& /*objectname*/ = "TpcSpaceChargeMatrixContainer" );

  /// calculate distortions by inverting stored matrices, and save relevant histograms
  void calculate_distortions();

  /// extrapolate distortions
  void extrapolate_distortions();

  /// save distortions
  void save_distortions(const std::string& /*filename*/ = "DistortionCorrections.root");

  //@}

  private:

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

  /// output distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_average;

  /// central membrane distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_cm;

};

#endif

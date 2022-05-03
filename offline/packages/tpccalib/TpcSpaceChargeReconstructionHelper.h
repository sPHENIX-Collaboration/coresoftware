#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
/**
 * \file TpcSpaceChargeReconstructionHelper.h
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>

// forward declarations
class TH3;
class TString;

class TpcSpaceChargeReconstructionHelper
{

  public:

  /// z extrapolation
  /**
   * interpolate between micromegas in the fully equiped sector
   */
  static void extrapolate_z( TH3* hin );

  /// get reference z range used for phi extrapolation normalization at a given radius
  static std::tuple<double, double> get_zref_range( double r );

  /// first phi extrapolation
  /**
   * copy the full z dependence of reference sector to all other sectors, separately for positive and negative z,
   * normalized by the measurement from provided micromegas, at the appropriate z
   */
  static void extrapolate_phi1( TH3* hin );

  /// second phi extrapolation
  /**
   * for each r, z and phi bin, linearly extrapolate between neighbor phi sector measurements
   */
  static void extrapolate_phi2( TH3* hin );

  /// separate positive and negative z histograms
  /**
   * split histograms in two, the first with negative z values only, the second with positive z values
   * this must be done before adding guarding bins around each axis, in order to prevent artifacts during calls to Interpolate
   * at the central membrane (z = 0)
   */
  static std::tuple<TH3*, TH3*> split( TH3* hin );

  /**
   * copy input histogram into output, with new name, while adding two "guarding bins" on
   * each axis, with identical content and error as the first and last bin of the original histogram
   * this is necessary for being able to call TH3->Interpolate() when using these histograms
   * to correct for the space charge distortions.
   * TODO: is this really necessary ? Possibly one could just use the bin content for the correction rather than using TH3->Interpolate,
   * in which case the "guarding bins" would be unnecessary. Should check if it leads to a significant deterioration of the momentum resolution
   */
  static TH3* copy_histogram( TH3* hin, const TString& name );

};

#endif

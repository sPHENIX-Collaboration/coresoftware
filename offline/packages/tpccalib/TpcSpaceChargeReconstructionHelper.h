#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
/**
 * \file TpcSpaceChargeReconstructionHelper.h
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>
#include <cstdint>

// forward declarations
class TH2;
class TH3;
class TString;

class TpcSpaceChargeReconstructionHelper
{
 public:
  /// create TPOT acceptance mask, using provided binning
  /**
   * this histogram contains 1 in cells that match the TPOT acceptance
   * only Micromegas in the central sector (TPC sectors 9 and 21) are considered
   */
  static void create_tpot_mask(TH3* /*source*/);

  /// first z extrapolation
  /**
   * interpolate along z between the two micromegas modules fully equiped sector
   * the mask is used to decide in which bins the extrapolation must be performed
   */
  static void extrapolate_z(TH3* /*source*/, const TH3* /*mask*/);

  /// second z extrapolation
  /**
   * interpolate along z between the readout plane (at which location the distortion is zero) and the outermost micromegas module
   * the mask is used to decide in which bins the extrapolation must be performed
   * the side is a bit mask to tell on which side of the histograms the pad-planes are
   */
  enum Side : uint8_t
  {
    Side_negative = 1 << 0,
    Side_positive = 1 << 1,
    Side_both = Side_negative | Side_positive
  };

  static void extrapolate_z2(TH3* /*source*/, const TH3* /*mask*/, uint8_t /* side */);

  /// first phi extrapolation
  /**
   * copy the full z dependence of reference sector to all other sectors
   * normalized by the measurement from provided micromegas, at the appropriate z
   * the mask is used to decide in which bins the extrapolation must be performed
   */
  static void extrapolate_phi1(TH3* /*source*/, const TH2* /*source_cm*/, const TH3* /*mask*/);

  /// second phi extrapolation
  /**
   * for each r, z and phi bin, linearly extrapolate between neighbor phi sector measurements
   */
  static void extrapolate_phi2(TH3* /*source*/, const TH3* /*mask*/);

  /// separate positive and negative z histograms
  /**
   * split histograms in two, the first with negative z values only, the second with positive z values
   * this must be done before adding guarding bins around each axis, in order to prevent artifacts during calls to Interpolate
   * at the central membrane (z = 0)
   */
  static std::tuple<TH3*, TH3*> split(const TH3* /*source*/);

  /**
   * copy input histogram into output, with new name, while adding two "guarding bins" on
   * each axis, with identical content and error as the first and last bin of the original histogram
   * this is necessary for being able to call TH3->Interpolate() when using these histograms
   * to correct for the space charge distortions.
   */
  static TH3* add_guarding_bins(const TH3* /*source*/, const TString& /*name*/);
};

#endif

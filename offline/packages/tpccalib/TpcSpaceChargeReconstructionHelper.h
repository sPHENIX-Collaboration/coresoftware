#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
/**
 * \file TpcSpaceChargeReconstructionHelper.h
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>
#include <cstdint>
#include <vector>

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

  /// shortcut to angular window, needed to define TPOT acceptance
  using range_t = std::pair<double, double>;

  /// list of angular windows
  using range_list_t = std::vector<range_t>;

  //!@name define phi range for each TPOT sector
  //@{
  static void set_phi_range_central( const range_t& );
  static void set_phi_range_east( const range_t& );
  static void set_phi_range_west( const range_t& );
  //@}

  //!@name define theta range in each TPOT sector.
  //@{
  static void set_theta_range_central( const range_list_t& range_list )
  { theta_range_central = range_list; }

  static void set_theta_range_east( const range_list_t& range_list )
  { theta_range_east = range_list; }

  static void set_theta_range_west( const range_list_t& range_list )
  { theta_range_west = range_list; }
  //@}

  private:

  ///@name detector geometry
  //@{
  /// phi range for central sector (TPC sectors 9 and 21)
  static range_t phi_range_central;

  /// phi range for east sector (TPC sectors 8 and 20)
  static range_t phi_range_east;

  /// phi range for west sector (TPC sectors 10 and 22)
  static range_t phi_range_west;

  /// list of theta angles for each micromeas in central sector with theta defined as atan2(z,r)
  static range_list_t theta_range_central;

  /// list of theta angles for each micromeas in central sector with theta defined as atan2(z,r)
  static range_list_t theta_range_east;

  /// list of theta angles for each micromeas in central sector with theta defined as atan2(z,r)
  static range_list_t theta_range_west;
  //@}

};

#endif

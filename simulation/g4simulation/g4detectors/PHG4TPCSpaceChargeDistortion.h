/// $Id: $

/*!
 * \file PHG4TPCSpaceChargeDistortion.h
 * \brief From TKH's SpaceChargeDistortion
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4TPCSPACECHARGEDISTORTION_H_
#define PHG4TPCSPACECHARGEDISTORTION_H_

class TH2D;

#include "PHG4TPCDistortion.h"

#include <string>

/*!
 * \brief PHG4TPCSpaceChargeDistortion
 */
///
///  Hello Space Charge DIstortion Fans:
///
///    This program is a quick & dirty bootstrap to allow reasonably accurate
///  space charge calculations to be made with minimal disturbance to our existing
///  TPC code base in sPHENIX.  Here is the concept of how it works:
///
///  To reasonable accuracy, we can model space charge distortions as having two
///  components, boith stemming dfrom the same original source:
///
///    1.  Loss of precision:
///        We assume that the precision loss is proportional to the DISTORTION
///        in a manner that is symmetric abouyt zero (after correction).
///    2.  Loss of accuracy:
///        Here we assume that we correctly poorly (wrong amount) and thereby
///        gain a SHIFT that is proportional the distortion itself.
///
///  To get a "best estimate" calculation, one uses an external analysis of space charge
///  distortions to determine the proportionality constants for each one of these.
///
///                                                                TKH
///                                                                5-19-2016
///
class PHG4TPCSpaceChargeDistortion : public PHG4TPCDistortion
{
public:
  PHG4TPCSpaceChargeDistortion(const std::string & distortion_map_file, int verbose = 0);

  virtual
  ~PHG4TPCSpaceChargeDistortion();

  //! radial distortion for a given truth location of the primary ionization
  double
  get_r_distortion(double r, double phi, double z) ;

  //! r*phi distortion for a given truth location of the primary ionization
  double
  get_rphi_distortion(double r, double phi, double z) ;

  //! z distortion for a given truth location of the primary ionization
  double
  get_z_distortion(double r, double phi, double z) {return 0;}

  TH2D *DRHIST()    {return rDistortion2;}
  TH2D *DRPHIHIST() {return rPhiDistortion2;}

  ///  NOTE:  ALICE claims that 20 cm deflections are corrected:
  ///         precision = 200 microns   (precisionFactor = 0.001)
  ///         accuracy = 1 mm           (accuracy factor = 0.005)
  ///  These defaults are hard-coded in the constructor, but the user
  ///  can change them using the routines below...
  void setPrecision(double p) {precisionFactor = p;}
  void setAccuracy (double a) { accuracyFactor = a;}

protected:

  TH2D *rDistortion2;
  TH2D *rPhiDistortion2;

  double precisionFactor;
  double accuracyFactor;

};

#endif /* PHG4TPCSPACECHARGEDISTORTION_H_ */

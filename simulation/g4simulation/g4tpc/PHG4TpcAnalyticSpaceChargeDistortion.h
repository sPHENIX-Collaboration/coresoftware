/// $Id: $

/*!
 * \file PHG4TPCAnalyticSpaceChargeDistortion.h
 * \brief Modified from TKH's SpaceChargeDistortion
 * \author Ross Corliss <ross.corliss@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4TPC_PHG4TPCANALYTICSPACECHARGEDISTORTION_H
#define G4TPC_PHG4TPCANALYTICSPACECHARGEDISTORTION_H

class TH2D;

#include "PHG4TpcDistortion.h"

#include <string>

/*!
 * \brief PHG4TpcAnalyticSpaceChargeDistortion
 */
///
///  Hello Space Charge Distortion Fans:
///
///    This program is a quick & dirty bootstrap to allow reasonably accurate
///  or completely inaccurate space charge distortions to be applied to hits
///  within the sPHENIX TPC simulation.  Note that this does not correct these
///  nor does it attempt to characterize the precision or accuracy of those
///  corrections.
///
///  We model the space charge as three TFormulas:
///    delr(r,phi,z), delphi(r,phi,z), and delz(r,phi,z)
///
///  These supply the signed shifts in r, phi, and z (time) positions at the
///  readout plane as a function of the r,phi, and z coordinates of the initial
///  particle position.
///
///  For realistic results, these can be pulled from fits to more detailed
///  calculations.
///
///                                                                -Ross
///                                                                15-Jan-2020
///
class TFormula;
class PHG4TpcAnalyticSpaceChargeDistortion : public PHG4TpcDistortion
{
 public:
  PHG4TpcAnalyticSpaceChargeDistortion( int verbose = 0);

  virtual ~PHG4TpcAnalyticSpaceChargeDistortion();

  //! radial distortion for a given truth location of the primary ionization
  double
    get_r_distortion(double r, double phi, double z);

  //! r*phi distortion for a given truth location of the primary ionization
  double
    get_rphi_distortion(double r, double phi, double z);

  //! z distortion for a given truth location of the primary ionization
  double
    get_z_distortion(double r, double phi, double z);

  
  void setFormulaR(TFormula form);
  void setFormulaPhi(TFormula form);
  void setFormulaZ(TFormula form);
  TFormula* getFormulaR();
  TFormula* getFormulaPhi();
  TFormula* getFormulaZ();

 protected:
  TFormula *delr;
  TFormula *delphi;
  TFormula *delz;

};

#endif /* G4TPC_PHG4TPCANALYTICSPACECHARGEDISTORTION_H */

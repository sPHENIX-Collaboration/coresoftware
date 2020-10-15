// $Id: $

/*!
 * \file PHG4TPCDistortion.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4TPC_PHG4TPCDISTORTION_H
#define G4TPC_PHG4TPCDISTORTION_H

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#include <TFile.h>
#include <TH3.h>
#include <TTree.h>
using namespace std;
#endif

/*!
 * \brief PHG4TpcDistortion is a virtual interface to apply distortion to a primary ionization in Tpc
 *
 * IO follow PHENIX units (cm)
 */


class TFile;
class PHG4TpcDistortion
{
 public:
  explicit PHG4TpcDistortion(int verbose = 0);

  virtual ~PHG4TpcDistortion();

  //! x distortion for a given truth location of the primary ionization
  double get_x_distortion(double x, double y, double z, int event_num);

  //! y distortion for a given truth location of the primary ionization
  double get_y_distortion(double x, double y, double z, int event_num);

  //! z distortion for a given truth location of the primary ionization
  double get_z_distortion(double x, double y, double z, int event_num);

  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void
  Verbosity(const int ival)
  {
    verbosity = ival;
  }

  //! Gets the verbosity of this module.
  virtual int
  Verbosity() const
  {
    return verbosity;
  }

  TH3F *hDXint;
  TH3F *hDYint;
  TH3F *hDZint;
  TH3F *TimehDZ;
  TH3F *TimehDX;
  TH3F *TimehDY;
  TTree *TimeTree;


 protected:
  //! The verbosity level. 0 means not verbose at all.
  int verbosity;

#if !defined(__CINT__) || defined(__CLING__)
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif /* G4TPC_PHG4TPCDISTORTION_H */

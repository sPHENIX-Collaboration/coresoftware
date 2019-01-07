// $Id: $

/*!
 * \file PHG4TPCDistortion.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4TPCDISTORTION_H_
#define PHG4TPCDISTORTION_H_

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

/*!
 * \brief PHG4TPCDistortion is a virtual interface to apply distortion to a primary ionization in TPC
 *
 * IO follow PHENIX units (cm)
 */
class PHG4TPCDistortion
{
 public:
  explicit PHG4TPCDistortion(int verbose = 0);

  virtual ~PHG4TPCDistortion();

  //! radial distortion for a given truth location of the primary ionization
  virtual double
  get_r_distortion(double r, double phi, double z) = 0;

  //! r*phi distortion for a given truth location of the primary ionization
  virtual double
  get_rphi_distortion(double r, double phi, double z) = 0;

  //! z distortion for a given truth location of the primary ionization
  virtual double
  get_z_distortion(double r, double phi, double z) = 0;

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

 protected:
  //! The verbosity level. 0 means not verbose at all.
  int verbosity;

#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif /* PHG4TPCDISTORTION_H_ */

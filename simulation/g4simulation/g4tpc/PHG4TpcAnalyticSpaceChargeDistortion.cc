// $Id: $

/*!
 * \file PHG4TpcAnalyticSpaceChargeDistortion.cc
 * \brief From TKH's SpaceChargeDistortion
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TpcAnalyticSpaceChargeDistortion.h"

#include <TAxis.h>            // for TAxis
#include <TObject.h>          // for TObject
#include <TFormula.h>

#include <gsl/gsl_randist.h>

#include <cassert>
#include <cstdlib>           // for exit, NULL
#include <iostream>

using namespace std;

PHG4TpcAnalyticSpaceChargeDistortion::PHG4TpcAnalyticSpaceChargeDistortion( int verbose)
  : PHG4TpcDistortion(verbose)
{
  //These formulae are the analytical forms of the change in that coordinate between
  //the point where the particle begins and the endcap.  They can be overwritten by individual calls
  //note that the canonical variable names, in order, are x,y,z, but they refer to r, phi, and z.

  float rmax=78.0;
  float rmin=20.0;
  
  delr=new TFormula("delr","[0]*cos([1]*(x+[2]))");
  //set the default secret numbers:
  delr->SetParameters(1.5,M_PI*2/(rmax-rmin),-rmin);
   delphi=new TFormula("delphi","0.0");
   delz=new TFormula("delz","0.0");
}

PHG4TpcAnalyticSpaceChargeDistortion::~PHG4TpcAnalyticSpaceChargeDistortion()
{
  if (delr)
    delete delr;
  if (delphi)
    delete delphi;
  if (delz)
    delete delz;
}

double
PHG4TpcAnalyticSpaceChargeDistortion::get_r_distortion(double r, double phi, double z)
{

  double del=delr->Eval(r,phi,z);
  if (verbosity > 0)
  {
    cout << "PHG4TpcAnalyticSpaceChargeDistortion::get_r_distortion - input"
         << " r = " << r
         << " phi = " << phi
         << " z = " << z << endl;
    cout << " R Distortion:" << del;
    cout << endl;
  }

  return del;
}

double
PHG4TpcAnalyticSpaceChargeDistortion::get_rphi_distortion(double r, double phi, double z)
{

  double del=r*delphi->Eval(r,phi,z);
  if (verbosity > 0)
  {
    cout << "PHG4TpcAnalyticSpaceChargeDistortion::get_r_distortion - input"
         << " r = " << r
         << " phi = " << phi
         << " z = " << z << endl;
    cout << " Rphi Distortion:" << del;
    cout << endl;
  }

  return del;
}

double
PHG4TpcAnalyticSpaceChargeDistortion::get_z_distortion(double r, double phi, double z)
{

  double del=delz->Eval(r,phi,z);
  if (verbosity > 0)
  {
    cout << "PHG4TpcAnalyticSpaceChargeDistortion::get_r_distortion - input"
         << " r = " << r
         << " phi = " << phi
         << " z = " << z << endl;
    cout << " Z Distortion:" << del;
    cout << endl;
  }

  return del;
}

// $Id: $

/*!
 * \file PHG4TpcDistortion.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>, modified by Henry Klest <henry.klest@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TpcDistortion.h"
#include "PHG4TpcElectronDrift.h"

#include <phool/PHRandomSeed.h>

#include <gsl/gsl_rng.h>         // for gsl_rng_alloc, gsl_rng_free, gsl_rng...
#include <cassert>
#include <cmath>                                       // for sqrt, fabs, NAN                                                                                                                               
#include <cstdlib>                                     // for exit                                                                                                                                          
#include <iostream>
#include <map>                                          // for _Rb_tree_cons...                                                                                                                             
#include <utility>                                      // for pair   
PHG4TpcDistortion::PHG4TpcDistortion(int verbose)
  : verbosity(verbose)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(RandomGenerator, seed);
  cout << "Inside PHG4TpcDistortion" << endl;
  TFile *StaticDistFile=new TFile("Summary_bX1508071_0_10_events.root.h_Charge_evt_0.real_B-1.5_E-400.0.ross_phi1_sphenix_phislice_lookup_r24xp24xz36.distortion_map.hist.root");//includes Trees of TH3Fs

  if(StaticDistFile->GetSize() == -1)
    {
      cout << "Distortion file could not be opened!" << endl;
    }
  //Open Static Space Charge Maps                                                                                                               
  hDXint=(TH3F*)StaticDistFile->Get("hIntDistortionX"); // Open TH3F files only once that contain distortions due to space charge                                                                     
  hDYint=(TH3F*)StaticDistFile->Get("hIntDistortionY");
  hDZint=(TH3F*)StaticDistFile->Get("hIntDistortionZ");

  TFile *TimeDistFile=new TFile("TimeOrderedDistortions.root");//includes Trees of TH3Fs                                                                                                                 
  if(TimeDistFile->GetSize() == -1)
    {
      cout << "TimeOrderedDistortion file could not be opened!" << endl;
    }
  TimeTree = (TTree*)TimeDistFile->Get("TimeDists");
  TimeTree->SetBranchAddress("hDistortionX",&TimehDX);
  TimeTree->SetBranchAddress("hDistortionY",&TimehDY);
  TimeTree->SetBranchAddress("hDistortionZ",&TimehDZ);

}
PHG4TpcDistortion::~PHG4TpcDistortion()
{
  gsl_rng_free(RandomGenerator);
}

double PHG4TpcDistortion::get_x_distortion(double x, double y, double z, int event_num)
{
  cout << "Inside PHG4TpcDistortion" << endl;
  double x_distortion;
  if(atan2(y,x) < 0)
    {
      x_distortion = hDXint->Interpolate(atan2(y,x)+2*M_PI,sqrt(x*x+y*y),z);
    }
  else
    {
      x_distortion = hDXint->Interpolate(atan2(y,x),sqrt(x*x+y*y),z);
    }
return x_distortion; 
}
double PHG4TpcDistortion::get_y_distortion(double x, double y, double z,int event_num)
{
double y_distortion;
  if(atan2(y,x) < 0)
    {
      y_distortion = hDYint->Interpolate(atan2(y,x)+2*M_PI,sqrt(x*x+y*y),z);
    }
  else
    {
      y_distortion = hDYint->Interpolate(atan2(y,x),sqrt(x*x+y*y),z);
    }
  return y_distortion; 
}
double PHG4TpcDistortion::get_z_distortion(double x, double y, double z, int event_num)
{
double z_distortion;
  if(atan2(y,x) < 0)
    {
      z_distortion = hDZint->Interpolate(atan2(y,x)+2*M_PI,sqrt(x*x+y*y),z);
    }
  else
    {
      z_distortion = hDZint->Interpolate(atan2(y,x),sqrt(x*x+y*y),z);
    }
  return z_distortion; 
}


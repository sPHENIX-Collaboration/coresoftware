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
PHG4TpcDistortion::PHG4TpcDistortion(int verbose,int event_num, bool do_time_ordered_distortion, bool do_static_distortion)
  : verbosity(verbose)
{
  //  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  // unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  // gsl_rng_set(RandomGenerator, seed);
 
 
  TFile *StaticDistFile=new TFile("$CALIBRATIONROOT/TPC/DistortionMaps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root");//includes Trees of TH3Fs
  if(do_static_distortion)
    {
  cout << "Using static TPC distortion map located at $CALIBRATIONROOT/TPC/DistortionMaps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root" << endl;
    }
  if(StaticDistFile->GetSize() == -1)
    {
      cout << "Distortion file could not be opened!" << endl;
    }
  //Open Static Space Charge Maps                                                                                                               
  hDXint=(TH3F*)StaticDistFile->Get("hIntDistortionX"); // Open TH3F files only once that contain distortions due to space charge                                                                     
  hDYint=(TH3F*)StaticDistFile->Get("hIntDistortionY");
  hDZint=(TH3F*)StaticDistFile->Get("hIntDistortionZ");
  
  TFile *TimeDistFile=new TFile("/gpfs/mnt/gpfs02/sphenix/user/klest/TimeOrderedDistortions.root");//includes Trees of TH3Fs                                                                
  // if(do_time_ordered_distortion)
  // {
  // cout << "Using Time-ordered TPC distortion map located at $CALIBRATIONROOT/TPC/DistortionMaps/TimeOrderedDistortions.root" << endl;
  // }                             
  if(TimeDistFile->GetSize() == -1)
    {
      cout << "TimeOrderedDistortion file could not be opened!" << endl;
    }
  TimehDX = new TH3F();
  TimehDY = new TH3F();
  TimehDZ = new TH3F();
  TimeTree = (TTree*)TimeDistFile->Get("TimeDists");
  TimeTree->SetBranchAddress("hIntDistortionX",&TimehDX);
  TimeTree->SetBranchAddress("hIntDistortionY",&TimehDY);
  TimeTree->SetBranchAddress("hIntDistortionZ",&TimehDZ);
  
}
PHG4TpcDistortion::~PHG4TpcDistortion()
{
  // gsl_rng_free(RandomGenerator);
}
void PHG4TpcDistortion::load_event(int event_num)
{
  TimeTree->GetEntry(event_num);
  return;
}

double PHG4TpcDistortion::get_x_distortion(double x, double y, double z,bool do_time_ordered_distortion,bool do_static_distortion)
{
  //cout << "TimehDX->GetZAxis()->FindBin(-50) is " << TimehDX->GetZaxis()->FindBin(-50) << " hDXInt bin -50 is " << hDXint->GetZaxis()->FindBin(-50) << endl;
  double phi = atan2(y,x);
  double r = sqrt(x*x+y*y);
  double x_distortion = 0;
  //cout << " do_static_distortion is " << do_static_distortion << " static x distortion is " << hDXint->Interpolate(atan2(y,x)+2*M_PI,r,z) << " do_time_ordered_distortion is " << do_time_ordered_distortion << " Time distortion is " << endl;//TimehDX->Interpolate(atan2(y,x)+2*M_PI,r,z) << endl;
  if(do_static_distortion && do_time_ordered_distortion) // brute forcing the boolean cases to be sure there's no undefined values sneaking in, there's definitely a better way to do this...
    {
      if(hDXint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  x_distortion = hDXint->Interpolate(phi+2*M_PI,r,z)+TimehDX->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  x_distortion = hDXint->Interpolate(phi,r,z)+TimehDX->Interpolate(phi,r,z);
	}
    }
  if(do_static_distortion && !do_time_ordered_distortion)
    {
      if(hDXint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  x_distortion = hDXint->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  x_distortion = hDXint->Interpolate(phi,r,z)+TimehDX->Interpolate(phi,r,z);
	}
    }
  if(!do_static_distortion && do_time_ordered_distortion)
    {
      if(phi < 0)
	{
	  x_distortion = TimehDX->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  x_distortion = TimehDX->Interpolate(phi,r,z)+TimehDX->Interpolate(phi,r,z);
	}
    }
      if(!do_time_ordered_distortion && !do_static_distortion)
	{x_distortion=0;}
return x_distortion; 
}
double PHG4TpcDistortion::get_y_distortion(double x, double y, double z,bool do_time_ordered_distortion,bool do_static_distortion)
{
  double phi = atan2(y,x);
  double r = sqrt(x*x+y*y); 
  double y_distortion = 0;
  if(do_static_distortion && do_time_ordered_distortion) // brute forcing the boolean cases to be sure there's no undefined values sneaking in, there's definitely a better way to do this... 
    {
      if(hDYint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  y_distortion = hDYint->Interpolate(phi+2*M_PI,r,z)+TimehDY->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  y_distortion = hDYint->Interpolate(phi,r,z)+TimehDY->Interpolate(phi,r,z)+TimehDX->Interpolate(phi,r,z);
	}
    }
  if(do_static_distortion && !do_time_ordered_distortion)
    {
      if(hDYint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  y_distortion = hDYint->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  y_distortion = hDYint->Interpolate(phi,r,z);
	}
    }
  if(!do_static_distortion && do_time_ordered_distortion)
    {
      if(phi < 0)
	{
	  y_distortion = TimehDY->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  y_distortion = TimehDY->Interpolate(phi,r,z);
	}
    }
  if(!do_time_ordered_distortion && !do_static_distortion)
    {y_distortion=0;}
  return y_distortion; 
}
double PHG4TpcDistortion::get_z_distortion(double x, double y, double z,bool do_time_ordered_distortion,bool do_static_distortion)
{
  double phi = atan2(y,x);
  double r = sqrt(x*x+y*y);
  double z_distortion = 0;
  if(do_static_distortion && do_time_ordered_distortion) // brute forcing the boolean cases to be sure there's no undefined values sneaking in, there's definitely a better way to do this... 
    {
      if(hDZint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  z_distortion = hDZint->Interpolate(phi+2*M_PI,r,z)+TimehDX->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  z_distortion = hDZint->Interpolate(phi,r,z)+TimehDX->Interpolate(phi,r,z);
	}
    }
  if(do_static_distortion && !do_time_ordered_distortion)
    {
      if(hDZint->GetZaxis()->FindBin(-50) == 0) // if z = -50 is in the underflow bin, map is only one-sided.
	{
	  z = fabs(z);
	}
      
      if(phi < 0)
	{
	  z_distortion = hDZint->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  z_distortion = hDZint->Interpolate(phi,r,z);
	}
    }
  if(!do_static_distortion && do_time_ordered_distortion)
    {
      if(phi < 0)
	{
	  z_distortion = TimehDZ->Interpolate(phi+2*M_PI,r,z);
	}
      else
	{
	  z_distortion = TimehDZ->Interpolate(phi,r,z);
	}
    }
  
if(!do_time_ordered_distortion && !do_static_distortion)
    {z_distortion=0;}
  return z_distortion; 
}


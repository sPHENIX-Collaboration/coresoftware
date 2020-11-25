// $Id: $

/*!
 * \file PHG4TpcDistortion.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>, Henry Klest <henry.klest@stonybrook.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TpcDistortion.h"

#include <TFile.h>
#include <TH3.h>
#include <TTree.h>

#include <cmath>                                       // for sqrt, fabs, NAN
#include <cstdlib>                                     // for exit
#include <iostream>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  //__________________________________________________________________________________
  double get_distortion(TH3* hstatic, TH3* htimeOrdered, double x, double y, double z)
  {

    double phi = std::atan2(y,x);
    if( phi < 0 ) phi += 2*M_PI;
    const double r = std::sqrt(square(x) + square(y));

    double x_distortion = 0;
    if( hstatic )
    {
      // if z = -50 is in the underflow bin, map is only one-sided.
      const auto zmap = ( hstatic->GetZaxis()->FindBin(-50) == 0) ? std::abs(z):z;
      x_distortion += hstatic->Interpolate( phi, r, zmap);
    }

    if( htimeOrdered )
    {
      // if z = -50 is in the underflow bin, map is only one-sided.
      const auto zmap = ( htimeOrdered->GetZaxis()->FindBin(-50) == 0) ? std::abs(z):z;
      x_distortion += htimeOrdered->Interpolate( phi, r, zmap);
    }

    return x_distortion;
  }

}  // namespace

//__________________________________________________________________________________________________________
PHG4TpcDistortion::PHG4TpcDistortion(bool do_time_ordered_distortion, bool do_static_distortion)
{
  if(do_static_distortion)
  {
    const TString filename( "$CALIBRATIONROOT/TPC/DistortionMaps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root");
    std::cout << "Using static TPC distortion map located at \"" << filename << "\"" << std::endl;

    m_staticTFile.reset( new TFile( filename ) );
    if( !m_staticTFile->IsOpen() )
    {
      std::cout << "Static distortion file could not be opened!" << std::endl;
      exit(1);
    }

    //Open Static Space Charge Maps
    hDXint = dynamic_cast<TH3*>(m_staticTFile->Get("hIntDistortionX"));
    hDYint = dynamic_cast<TH3*>(m_staticTFile->Get("hIntDistortionY"));
    hDZint = dynamic_cast<TH3*>(m_staticTFile->Get("hIntDistortionZ"));
  }

  if(do_time_ordered_distortion)
  {
    const TString filename("/gpfs/mnt/gpfs02/sphenix/user/klest/TimeOrderedDistortions.root");
    m_timeOrderedTFile.reset( new TFile( filename ) );
    if( !m_timeOrderedTFile->IsOpen() )
    {
      std::cout << "TimeOrdered distortion file could not be opened!" << std::endl;
      exit(1);
    }

    // create histograms
    TimehDX = new TH3F();
    TimehDY = new TH3F();
    TimehDZ = new TH3F();

    TimeTree = static_cast<TTree*>( m_timeOrderedTFile->Get("TimeDists") );
    TimeTree->SetBranchAddress("hIntDistortionX",&TimehDX);
    TimeTree->SetBranchAddress("hIntDistortionY",&TimehDY);
    TimeTree->SetBranchAddress("hIntDistortionZ",&TimehDZ);
  }

}

//__________________________________________________________________________________________________________
void PHG4TpcDistortion::load_event(int event_num)
{

  if( TimeTree )
  {
    int nentries = TimeTree->GetEntries();
    if(event_num > nentries) event_num = event_num % nentries;
    if(event_num%nentries == 0 && event_num != 0)
    {
      std::cout << "Distortion map sequence repeating as of event number " << event_num  << std::endl;
    }
    TimeTree->GetEntry(event_num);
  }

  return;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_x_distortion(double x, double y, double z)
{ return get_distortion(hDXint, TimehDX, x, y, z ); }

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_y_distortion(double x, double y, double z)
{ return get_distortion(hDYint, TimehDY, x, y, z ); }

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion(double x, double y, double z)
{ return get_distortion(hDZint, TimehDZ, x, y, z ); }

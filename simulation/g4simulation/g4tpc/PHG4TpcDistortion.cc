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

#include <cmath>    // for sqrt, fabs, NAN
#include <cstdlib>  // for exit
#include <iostream>

namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  // check boundaries in axis
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries( const TAxis* axis, double value )
  {
    const auto bin = axis->FindBin( value );
    return( bin >= 2 && bin < axis->GetNbins() );
  }
  
  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries( const TH3* h, double r, double phi, double z )
  {
    return check_boundaries( h->GetXaxis(), r ) 
      && check_boundaries( h->GetYaxis(), phi ) 
      && check_boundaries( h->GetZaxis(), z );
  }
  
}  // namespace

//__________________________________________________________________________________________________________
void PHG4TpcDistortion::Init()
{
  if (m_do_static_distortions)
  {
    std::cout << "PHG4TpcDistortion::Init - m_static_distortion_filename: " << m_static_distortion_filename << std::endl;
    m_static_tfile.reset(new TFile(m_static_distortion_filename.c_str()));
    if (!m_static_tfile->IsOpen())
    {
      std::cout << "Static distortion file could not be opened!" << std::endl;
      exit(1);
    }

    //Open Static Space Charge Maps
    hDRint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionR_negz"));
    hDRint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionR_posz"));
    hDPint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionP_negz"));
    hDPint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionP_posz"));
    hDZint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionZ_negz"));
    hDZint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionZ_posz"));
  }

  if (m_do_time_ordered_distortions)
  {
    std::cout << "PHG4TpcDistortion::Init - m_time_ordered_distortion_filename: " << m_time_ordered_distortion_filename << std::endl;
    m_time_ordered_tfile.reset(new TFile(m_time_ordered_distortion_filename.c_str()));
    if (!m_time_ordered_tfile->IsOpen())
    {
      std::cout << "TimeOrdered distortion file could not be opened!" << std::endl;
      exit(1);
    }

    // create histograms
    TimehDR[0] = new TH3F();
    TimehDR[1] = new TH3F();
    TimehDP[0] = new TH3F();
    TimehDP[1] = new TH3F();
    TimehDZ[0] = new TH3F();
    TimehDZ[1] = new TH3F();

    TimeTree = static_cast<TTree*>(m_time_ordered_tfile->Get("TimeDists"));
    TimeTree->SetBranchAddress("hIntDistortionR_negz", &(TimehDR[0]));
    TimeTree->SetBranchAddress("hIntDistortionR_posz", &(TimehDR[1]));
    TimeTree->SetBranchAddress("hIntDistortionP_negz", &(TimehDP[0]));
    TimeTree->SetBranchAddress("hIntDistortionP_posz", &(TimehDP[1]));
    TimeTree->SetBranchAddress("hIntDistortionZ_negz", &(TimehDZ[0]));
    TimeTree->SetBranchAddress("hIntDistortionZ_posz", &(TimehDZ[1]));
  }
}

//__________________________________________________________________________________________________________
void PHG4TpcDistortion::load_event(int event_num)
{
  if (TimeTree)
  {
    int nentries = TimeTree->GetEntries();
    if (event_num > nentries) event_num = event_num % nentries;
    if (event_num % nentries == 0 && event_num != 0)
    {
      std::cout << "Distortion map sequence repeating as of event number " << event_num << std::endl;
    }
    TimeTree->GetEntry(event_num);
  }

  return;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_x_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  //get components
  double dr = get_distortion('r', r, phi, z);
  double dphi = get_distortion('p', r, phi, z);

  //rotate into cartesian based on local r phi:
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double dx = dr * cosphi - dphi * sinphi;
  return dx;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_y_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  //get components
  double dr = get_distortion('r', r, phi, z);
  double dphi = get_distortion('p', r, phi, z);

  //rotate into cartesian based on local r phi:
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double dy = dphi * cosphi + dr * sinphi;
  return dy;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  //get components
  double dz = get_distortion('z', r, phi, z);

  return dz;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_r_distortion(double r, double phi, double z) const
{
  return get_distortion('r', r, phi, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_rphi_distortion(double r, double phi, double z) const
{
  return get_distortion('p', r, phi, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion(double r, double phi, double z) const
{
  return get_distortion('z', r, phi, z);
}

double PHG4TpcDistortion::get_distortion(char axis, double r, double phi, double z) const
{
  if (phi < 0) phi += 2 * M_PI;
  const int zpart = (z > 0 ? 1 : 0);  //z<0 corresponds to the negative side, which is element 0.

  TH3* hdistortion = nullptr;

  if (axis != 'r' && axis != 'p' && axis != 'z')
  {
    std::cout << "Distortion Requested along axis " << axis << " which is invalid.  Exiting.\n"
              << std::endl;
    exit(1);
  }

  double _distortion = 0.;

  //select the appropriate histogram:
  if (m_do_static_distortions)
  {
    if (axis == 'r')
    {
      hdistortion = hDRint[zpart];
    }
    else if (axis == 'p')
    {
      hdistortion = hDPint[zpart];
    }
    else if (axis == 'z')
    {
      hdistortion = hDZint[zpart];
    }
    if (hdistortion)
    {
      if( check_boundaries( hdistortion, phi, r, z ) )
      { _distortion += hdistortion->Interpolate(phi, r, z); }
    }
    else
    {
      std::cout << "Static Distortion Requested along axis " << axis << ", but distortion map does not exist.  Exiting.\n"
                << std::endl;
      exit(1);
    }
  }

  if (m_do_time_ordered_distortions)
  {
    if (axis == 'r')
    {
      hdistortion = TimehDR[zpart];
    }
    else if (axis == 'p')
    {
      hdistortion = TimehDP[zpart];
    }
    else if (axis == 'z')
    {
      hdistortion = TimehDZ[zpart];
    }
    if (hdistortion)
    {
      if( check_boundaries( hdistortion, phi, r, z ) )
      { _distortion += hdistortion->Interpolate(phi, r, z); }
    }
    else
    {
      std::cout << "Time Series Distortion Requested along axis " << axis << ", but distortion map does not exist.  Exiting.\n"
                << std::endl;
      exit(1);
    }
  }

  return _distortion;
}

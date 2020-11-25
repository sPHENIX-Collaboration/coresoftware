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

  //__________________________________________________________________________________
  double get_distortion(TH3* hstatic, TH3* htimeOrdered, double x, double y, double z)
  {
    double phi = std::atan2(y, x);
    if (phi < 0) phi += 2 * M_PI;
    const double r = std::sqrt(square(x) + square(y));

    double x_distortion = 0;
    if (hstatic)
    {
      // if z = -50 is in the underflow bin, map is only one-sided.
      const auto zmap = (hstatic->GetZaxis()->FindBin(-50) == 0) ? std::abs(z) : z;
      x_distortion += hstatic->Interpolate(phi, r, zmap);
    }

    if (htimeOrdered)
    {
      // if z = -50 is in the underflow bin, map is only one-sided.
      const auto zmap = (htimeOrdered->GetZaxis()->FindBin(-50) == 0) ? std::abs(z) : z;
      x_distortion += htimeOrdered->Interpolate(phi, r, zmap);
    }

    return x_distortion;
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
    hDXint = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionX"));
    hDYint = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionY"));
    hDZint = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionZ"));
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
    TimehDX = new TH3F();
    TimehDY = new TH3F();
    TimehDZ = new TH3F();

    TimeTree = static_cast<TTree*>(m_time_ordered_tfile->Get("TimeDists"));
    TimeTree->SetBranchAddress("hIntDistortionX", &TimehDX);
    TimeTree->SetBranchAddress("hIntDistortionY", &TimehDY);
    TimeTree->SetBranchAddress("hIntDistortionZ", &TimehDZ);
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
double PHG4TpcDistortion::get_x_distortion(double x, double y, double z)
{
  return get_distortion(hDXint, TimehDX, x, y, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_y_distortion(double x, double y, double z)
{
  return get_distortion(hDYint, TimehDY, x, y, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion(double x, double y, double z)
{
  return get_distortion(hDZint, TimehDZ, x, y, z);
}

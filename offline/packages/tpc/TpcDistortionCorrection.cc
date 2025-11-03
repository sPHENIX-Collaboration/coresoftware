/*!
 * \file TpcDistortionCorrection.cc
 * \brief applies provided distortion corrections to a cluster and returns corrected position
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDistortionCorrection.h"
#include "TpcDistortionCorrectionContainer.h"

#include <TH1.h>
#include <cmath>

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
  inline bool check_boundaries(const TAxis* axis, double value)
  {
    const auto bin = axis->FindBin(value);
    return (bin >= 2 && bin < axis->GetNbins());
  }

  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TH1* h, double r, double phi, double z)
  {
    return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi) && check_boundaries(h->GetZaxis(), z);
  }

  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TH1* h, double r, double phi)
  {
    return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi);
  }

}  // namespace

//________________________________________________________
Acts::Vector3 TpcDistortionCorrection::get_corrected_position(const Acts::Vector3& source, const TpcDistortionCorrectionContainer* dcc, unsigned int mask) const
{
  // get cluster radius, phi and z
  const auto r = std::sqrt(square(source.x()) + square(source.y()));
  auto phi = std::atan2(source.y(), source.x());
  if (phi < 0)
  {
    phi += 2 * M_PI;
  }

  const auto z = source.z();
  const int index = z > 0 ? 1 : 0;

  // apply corrections
  auto phi_new = phi;
  auto r_new = r;
  auto z_new = z;

  // if the phi correction hist units are cm, we must divide by r to get the dPhi in radians
  auto divisor = r;

  if (dcc->m_phi_hist_in_radians)
  {
    // if the phi correction hist units are radians, we must not divide by r.
    divisor = 1.0;
  }

  //set our default corrections to be zero:
  //first inherit the same type
  auto dphi=phi;
  auto dr=r;
  auto dz=z;
  //then set them to zero:
  dphi=0;
  dr=0;
  dz=0;
  
  //get the corrections from the histograms
  if (dcc->m_dimensions == 3)
  {
    if (dcc->m_hDPint[index] && (mask & COORD_PHI) && check_boundaries(dcc->m_hDPint[index], phi, r, z))
    {
      dphi=dcc->m_hDPint[index]->Interpolate(phi, r, z) / divisor;
    }
    if (dcc->m_hDRint[index] && (mask & COORD_R) && check_boundaries(dcc->m_hDRint[index], phi, r, z))
    {
      dr=dcc->m_hDRint[index]->Interpolate(phi, r, z);
    }
    if (dcc->m_hDZint[index] && (mask & COORD_Z) && check_boundaries(dcc->m_hDZint[index], phi, r, z))
    {
      dz=dcc->m_hDZint[index]->Interpolate(phi, r, z);
    }
  }
  else if (dcc->m_dimensions == 2)
  {
    double zterm = 1.0;

    if (dcc->m_interpolate_z){
      zterm=(1. - std::abs(z) / 102.605);
    }
    if (dcc->m_hDPint[index] && (mask & COORD_PHI) && check_boundaries(dcc->m_hDPint[index], phi, r))
    {
      dphi=dcc->m_hDPint[index]->Interpolate(phi, r) * zterm / divisor;
    }
    if (dcc->m_hDRint[index] && (mask & COORD_R) && check_boundaries(dcc->m_hDRint[index], phi, r))
    {
      dr=dcc->m_hDRint[index]->Interpolate(phi, r) * zterm;
    }
    if (dcc->m_hDZint[index] && (mask & COORD_Z) && check_boundaries(dcc->m_hDZint[index], phi, r))
    {
      dz=dcc->m_hDZint[index]->Interpolate(phi, r) * zterm;
    }
    
  }

//if we are scaling, apply the scale factor to each correction
if(dcc->m_use_scalefactor)
  {
    dphi *= dcc->m_scalefactor;
    dr *= dcc->m_scalefactor;
    dz *= dcc->m_scalefactor;
  }

  phi_new=phi-dphi;
  r_new=r-dr;
  z_new=z-dz;

  // update cluster
  const auto x_new = r_new * std::cos(phi_new);
  const auto y_new = r_new * std::sin(phi_new);

  return {x_new, y_new, z_new};
}

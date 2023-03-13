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
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // check boundaries in axis
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries( const TAxis* axis, double value )
  {
    const auto bin = axis->FindBin( value );
    return( bin >= 2 && bin < axis->GetNbins() );
  }
  
  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries( const TH1* h, double r, double phi, double z )
  {
    return check_boundaries( h->GetXaxis(), r ) 
      && check_boundaries( h->GetYaxis(), phi ) 
      && check_boundaries( h->GetZaxis(), z );
  }
 
  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries( const TH1* h, double r, double phi )
  { return check_boundaries( h->GetXaxis(), r ) && check_boundaries( h->GetYaxis(), phi ); }

}

//________________________________________________________
Acts::Vector3 TpcDistortionCorrection::get_corrected_position( const Acts::Vector3& source, const TpcDistortionCorrectionContainer* dcc, unsigned int mask) const
{
  // get cluster radius, phi and z
  const auto r = std::sqrt( square( source.x() ) + square( source.y() ) );
  auto phi = std::atan2( source.y(), source.x() );
  if( phi < 0 ) phi += 2*M_PI;

  const auto z = source.z();
  const int index = z > 0 ? 1:0;

  // apply corrections
  auto phi_new=phi;
  auto r_new=r;
  auto z_new=z;
  if (dcc->dimensions==3)
  {
    if (dcc->m_hDPint[index] && (mask&COORD_PHI) && check_boundaries( dcc->m_hDPint[index],phi,r,z)) phi_new = phi - dcc->m_hDPint[index]->Interpolate(phi,r,z)/r;
    if (dcc->m_hDRint[index] && (mask&COORD_R) && check_boundaries( dcc->m_hDRint[index],phi,r,z)) r_new = r - dcc->m_hDRint[index]->Interpolate(phi,r,z);
    if (dcc->m_hDZint[index] && (mask&COORD_Z) && check_boundaries( dcc->m_hDZint[index],phi,r,z)) z_new = z - dcc->m_hDZint[index]->Interpolate(phi,r,z);
  }
  else if (dcc->dimensions==2){
    const double zterm = (1.- std::abs(z)/105.5);
    if (dcc->m_hDPint[index] && (mask&COORD_PHI) && check_boundaries( dcc->m_hDPint[index],phi,r)) phi_new = phi - dcc->m_hDPint[index]->Interpolate(phi,r)*zterm/r;
    if (dcc->m_hDRint[index] && (mask&COORD_R) && check_boundaries( dcc->m_hDRint[index],phi,r)) r_new = r - dcc->m_hDRint[index]->Interpolate(phi,r)*zterm;
    if (dcc->m_hDZint[index] && (mask&COORD_Z) && check_boundaries( dcc->m_hDZint[index],phi,r)) z_new = z - dcc->m_hDZint[index]->Interpolate(phi,r)*zterm;
  }
  
  // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  return {x_new, y_new, z_new};
}

/*!
 * \file TpcDistortionCorrection.cc
 * \brief applies provided distortion corrections to a cluster and returns corrected position
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDistortionCorrection.h"
#include "TpcDistortionCorrectionObject.h"

#include <TH3.h>
#include <cmath>

namespace
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
}

//________________________________________________________
TVector3 TpcDistortionCorrection::get_corrected_position( const TVector3& source, const TpcDistortionCorrectionObject* dco) const
{
  // get cluster radius, phi and z
  const auto r = std::sqrt( square( source.x() ) + square( source.y() ) );
  auto phi = std::atan2( source.y(), source.x() );
  if( phi < 0 ) phi += 2*M_PI;

  const auto z = source.z();
  const int index = z > 0 ? 1:0;

  // apply corrections
  const auto phi_new = phi - dco->m_hDPint[index]->Interpolate(phi,r,z)/r;
  const auto r_new = r - dco->m_hDRint[index]->Interpolate(phi,r,z);
  const auto z_new = z - dco->m_hDZint[index]->Interpolate(phi,r,z);
  
  // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  return {x_new, y_new, z_new};
}

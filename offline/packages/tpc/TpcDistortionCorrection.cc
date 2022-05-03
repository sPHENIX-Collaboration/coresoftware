/*!
 * \file TpcDistortionCorrection.cc
 * \brief applies provided distortion corrections to a cluster and returns corrected position
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDistortionCorrection.h"
#include "TpcDistortionCorrectionContainer.h"

#include <TH3.h>
#include <cmath>

namespace
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
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
  if (dcc->dimensions==3){
    phi_new = (dcc->m_hDPint[index] && (mask&COORD_PHI)) ? phi - dcc->m_hDPint[index]->Interpolate(phi,r,z)/r : phi;
    r_new = (dcc->m_hDRint[index] && (mask&COORD_R)) ? r - dcc->m_hDRint[index]->Interpolate(phi,r,z) : r;
    z_new = (dcc->m_hDZint[index] && (mask&COORD_Z)) ? z - dcc->m_hDZint[index]->Interpolate(phi,r,z) : z;
  }
  else if (dcc->dimensions==2){
    const auto bin=dcc->m_hDPint[index]->FindBin(phi,r,z);
    const auto zterm=(1-z/105.5); //where to get the official z position limits?
    phi_new = (dcc->m_hDPint[index] && (mask&COORD_PHI)) ? phi - dcc->m_hDPint[index]->GetBinContent(bin)/r*zterm : phi;
    r_new = (dcc->m_hDRint[index] && (mask&COORD_R)) ? r - dcc->m_hDRint[index]->GetBinContent(bin)*zterm : r;
    z_new = (dcc->m_hDZint[index] && (mask&COORD_Z)) ? z - dcc->m_hDZint[index]->GetBinContent(bin)*zterm : z;
  }
  else {
    //crash?
  }
  
    // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  return {x_new, y_new, z_new};
}

/**
 * @file trackbase/CMFlashDifferencev1.cc
 * @author Tony Frawley
 * @date January 2022
 * @brief Implementation of CMFlashDifferencev1
 */
#include "CMFlashDifferencev1.h"

#include <cmath>
#include <utility>          // for swap

void CMFlashDifferencev1::identify(std::ostream& os) const
{
  os << "---CMFlashDifferencev1--------------------" << std::endl;

  os << "nclusters: " << m_nclusters << std::dec << std::endl;

  os << " truth Phi =  " << m_Phi[0];
  os << ",  reco Phi = " << m_Phi[1] << ") rad";

  os << " truth R =  " << m_R[0];
  os << ",  reco R = " << m_R[1] << ") cm";

  os << " truth Z =  " << m_Z[0];
  os << ",  reco Z = " << m_Z[1] << ") cm";

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int CMFlashDifferencev1::isValid() const
{
  if (m_nclusters == UINT_MAX) return 0;

  if(std::isnan(getTruthPhi())) return 0;
  if(std::isnan(getTruthR())) return 0;
  if(std::isnan(getTruthZ())) return 0;

  if(std::isnan(getRecoPhi())) return 0;
  if(std::isnan(getRecoR())) return 0;
  if(std::isnan(getRecoZ())) return 0;


  return 1;
}

void CMFlashDifferencev1::CopyFrom( const CMFlashDifference& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  CMFlashDifference::CopyFrom( source );

  setNclusters( source.getNclusters() );

  setTruthPhi( source.getTruthPhi() );
  setTruthR( source.getTruthR() );
  setTruthZ( source.getTruthZ() );

  setRecoPhi( source.getRecoPhi() );
  setRecoR( source.getRecoR() );
  setRecoZ( source.getRecoZ() );

}


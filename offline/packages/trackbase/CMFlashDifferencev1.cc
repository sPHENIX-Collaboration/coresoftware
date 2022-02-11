/**
 * @file trackbase/CMFlashDifferencev1.cc
 * @author Tony Frawley
 * @date January 2022
 * @brief Implementation of CMFlashDifferencev1
 */
#include "CMFlashDifferencev1.h"

#include <cmath>
#include <utility>          // for swap

namespace
{

}

CMFlashDifferencev1::CMFlashDifferencev1()
  : m_key(UINT_MAX)
{

  for (int i = 0; i < 2; ++i) m_pos_truth[i] = NAN;
  for (int i = 0; i < 2; ++i) m_pos_reco[i] = NAN;

 }

void CMFlashDifferencev1::identify(std::ostream& os) const
{
  os << "---CMFlashDifferencev1--------------------" << std::endl;
  os << "key: " << getKey() << std::dec << std::endl;

  os << " truth (x,y) =  (" << m_pos_truth[0];
  os << ", " << m_pos_truth[1] << ") cm";

  os << " reco (x,y) =  (" << m_pos_reco[0];
  os << ", " << m_pos_reco[1] << ") cm";


  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int CMFlashDifferencev1::isValid() const
{
  if (m_key == UINT_MAX) return 0;

  if(std::isnan(getTruthX())) return 0;
  if(std::isnan(getTruthY())) return 0;
  if(std::isnan(getRecoX())) return 0;
  if(std::isnan(getRecoY())) return 0;

  return 1;
}

void CMFlashDifferencev1::CopyFrom( const CMFlashDifference& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  CMFlashDifference::CopyFrom( source );

  setKey( source.getKey() );
  setTruthX( source.getTruthX() );
  setTruthY( source.getTruthY() );
  setRecoX( source.getRecoX() );
  setRecoY( source.getRecoY() );

}


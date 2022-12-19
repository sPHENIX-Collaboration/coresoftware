/**
 * @file trackbase/RawHitSetv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of RawHitSetv1
 */
#include "RawHitSetv1.h"
#include "RawHit.h"

#include <cstdlib>     // for exit
#include <iostream>
#include <type_traits>  // for __decay_and_strip<>::__type

void RawHitSetv1::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;
  for( const auto& hit:m_hits )
  { delete hit; }

  m_hits.clear();
  return;
}

void RawHitSetv1::identify(std::ostream& os) const
{
  const unsigned int layer = TrkrDefs::getLayer(m_hitSetKey);
  const unsigned int trkrid =  TrkrDefs::getTrkrId(m_hitSetKey);    
  os 
    << "RawHitSetv1: "   
    << "       hitsetkey " << getHitSetKey()
    << " TrkrId " << trkrid 
    << " layer " << layer
    << " nhits: " << m_hits.size() 
    << std::endl;

}

void RawHitSetv1::addHit(RawHit* hit)
{
  m_hits.push_back(hit);
  return;
}
//void RawHitSetv1::addTpcHit(unsigned short phibin,RawHit* hit)
///{
//  m_tpchits[phibin].push_back(hit);
// return;
//}
void RawHitSetv1::setTpcPhiBins(unsigned short phibins)
{
  m_tpchits.resize(phibins);
  return;
}
/*RawHit* 
RawHitSetv1::getHit(const TrkrDefs::hitkey key) const
{
  RawHitSetv1::ConstIterator it = m_hits.find(key);
  
  if (it != m_hits.end()) return it->second;
  else return nullptr;
}
*/
RawHitSetv1::ConstRange RawHitSetv1::getHits() const
{
  return std::make_pair(m_hits.cbegin(), m_hits.cend());
}

//RawHitSetv1::ConstRange RawHitSetv1::getTpcHits(unsigned short phibin) const
//{
//  return std::make_pair(m_tpchits[phibin].cbegin(), m_tpchits[phibin].cend());
//}

/**
 * @file trackbase/TrkrHitSetv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitSetv1
 */
#include "TrkrHitSetv1.h"
#include "TrkrHit.h"

#include <cstdlib>     // for exit
#include <iostream>
#include <type_traits>  // for __decay_and_strip<>::__type

void TrkrHitSetv1::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  for(auto&& [key, hit] : m_hits)
    { delete hit; }

  m_hits.clear();
  
}

void TrkrHitSetv1::identify(std::ostream& os) const
{
  const unsigned int layer = TrkrDefs::getLayer(m_hitSetKey);
  const unsigned int trkrid =  TrkrDefs::getTrkrId(m_hitSetKey);    
  os 
    << "TrkrHitSetv1: "   
    << "       hitsetkey " << getHitSetKey()
    << " TrkrId " << trkrid 
    << " layer " << layer
    << " nhits: " << m_hits.size() 
    << std::endl;

  for( const auto& entry : m_hits )
  {
    std::cout << " hitkey " << entry.first << std::endl;
    (entry.second)->identify(os);
  }
}

void TrkrHitSetv1::removeHit(TrkrDefs::hitkey key)
{
  const auto it = m_hits.find(key);
  if (it != m_hits.end())
  {
    delete it->second;
    m_hits.erase(it);
  } else {
    identify();
    std::cout << "TrkrHitSetv1::removeHit: deleting a nonexist key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}

TrkrHitSetv1::ConstIterator
TrkrHitSetv1::addHitSpecificKey(const TrkrDefs::hitkey key, TrkrHit* hit)
{
  const auto ret = m_hits.insert(std::make_pair(key, hit));
  if ( !ret.second )
  {
    std::cout << "TrkrHitSetv1::AddHitSpecificKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  } else {
    return ret.first;
  }
}

TrkrHit* 
TrkrHitSetv1::getHit(const TrkrDefs::hitkey key) const
{
  TrkrHitSetv1::ConstIterator it = m_hits.find(key);
  
  if (it != m_hits.end()) return it->second;
  else return nullptr;
}

TrkrHitSetv1::ConstRange
TrkrHitSetv1::getHits() const
{
  return std::make_pair(m_hits.cbegin(), m_hits.cend());
}

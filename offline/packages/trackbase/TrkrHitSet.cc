/**
 * @file trackbase/TrkrHitSet.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitSet
 */
#include "TrkrHitSet.h"

#include "TrkrHit.h"

#include <cstdlib>     // for exit
#include <iostream>
#include <type_traits>  // for __decay_and_strip<>::__type


TrkrHitSet::TrkrHitSet()
  : m_hitSetKey(TrkrDefs::HITSETKEYMAX)
  , m_hits()
{
}

TrkrHitSet::~TrkrHitSet()
{
  // frees the memory associated with the pointers gto the TrkrHit objects
  Reset();
}

void TrkrHitSet::print() const
{
  identify(std::cout);
}

void TrkrHitSet::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;
  
  while ( m_hits.begin() != m_hits.end() )
  {
    delete (m_hits.begin())->second;
    m_hits.erase(m_hits.begin());
  }
  
  return;
}

void TrkrHitSet::identify(std::ostream& os) const
{
  int layer = TrkrDefs::getLayer(getHitSetKey());
  int trkrid =  TrkrDefs::getTrkrId(getHitSetKey());    
  os << "TrkrHitSet: "   << "       hitsetkey " << getHitSetKey() << " TrkrId " << trkrid << " layer " << layer << " nhits: " << m_hits.size() << std::endl;

  for ( auto& entry : m_hits )
  {
    std::cout << " hitkey " << entry.first << std::endl;
    (entry.second)->identify(os);
  }
}


void TrkrHitSet::removeHit(TrkrDefs::hitkey key)
{

  TrkrHitSet::ConstIterator it = m_hits.find(key);

  if (it != m_hits.end())
  {
    delete it->second;
    m_hits.erase(key);
  }
  else
  {
    identify();
    std::cout << "TrkrHitSet::removeHit: deleting a nonexist key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}



TrkrHitSet::ConstIterator
TrkrHitSet::addHitSpecificKey(const TrkrDefs::hitkey key, TrkrHit* hit)
{
  auto ret = m_hits.insert(std::make_pair(key, hit));

  if ( !ret.second )
  {
    std::cout << "TrkrHitSet::AddHitSpecificKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

TrkrHit*
TrkrHitSet::getHit(const TrkrDefs::hitkey key)
{
  TrkrHitSet::ConstIterator it = m_hits.find(key);
  
  if (it != m_hits.end())
  {
    return it->second;
  }
  
  return nullptr;
}

TrkrHitSet::ConstRange
TrkrHitSet::getHits()
{
  return std::make_pair(m_hits.begin(), m_hits.end());
}

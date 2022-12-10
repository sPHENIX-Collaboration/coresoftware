/**
 * @file trackbase/RawHitSetContainerv1.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation for RawHitSetContainerv1
 */
#include "RawHitSetContainerv1.h"

#include "TrkrDefs.h"
#include "RawHitSetv1.h"

#include <cstdlib>

void RawHitSetContainerv1::Reset()
{
  for( const auto& pair:m_hitmap )
  { delete pair.second; }
  
  m_hitmap.clear();
}

void RawHitSetContainerv1::identify(std::ostream& os) const
{
  ConstIterator iter;
  os << "Number of hits: " << size() << std::endl;
  for( const auto& pair:m_hitmap )
  {
    int layer = TrkrDefs::getLayer(pair.first);
    os << "hitsetkey " << pair.first << " layer " <<layer << std::endl;
    pair.second->identify();
  }
  return;
}

RawHitSetContainerv1::ConstIterator
RawHitSetContainerv1::addHitSet(RawHitSet* newhit)
{ return addHitSetSpecifyKey(newhit->getHitSetKey(), newhit); }

RawHitSetContainerv1::ConstIterator
RawHitSetContainerv1::addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, RawHitSet* newhit)
{
  const auto ret = m_hitmap.insert(std::make_pair(key, newhit));
  if ( !ret.second )
  {
    std::cout << "RawHitSetContainerv1::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  } else {
    return ret.first;
  }
}

void RawHitSetContainerv1::removeHitSet(TrkrDefs::hitsetkey key)
{ m_hitmap.erase(key); }

void RawHitSetContainerv1::removeHitSet(RawHitSet *hitset)
{ removeHitSet( hitset->getHitSetKey() ); }

RawHitSetContainerv1::ConstRange
RawHitSetContainerv1::getHitSets(const TrkrDefs::TrkrId trackerid) const
{
  const TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid);
  const TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid);
  return std::make_pair( m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi) );
}

RawHitSetContainerv1::ConstRange
RawHitSetContainerv1::getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
{
  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid, layer);
  return std::make_pair( m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi) );
}

RawHitSetContainerv1::ConstRange
RawHitSetContainerv1::getHitSets() const
{ return std::make_pair(m_hitmap.cbegin(), m_hitmap.cend()); }

RawHitSetContainerv1::Iterator
RawHitSetContainerv1::findOrAddHitSet(TrkrDefs::hitsetkey key)
{
  auto it = m_hitmap.lower_bound( key );
  if( it == m_hitmap.end() || (key < it->first ) )
  {
    it = m_hitmap.insert(it, std::make_pair(key, new RawHitSetv1));
    it->second->setHitSetKey( key );
  }
  return it;
}

RawHitSet*
RawHitSetContainerv1::findHitSet(TrkrDefs::hitsetkey key)
{
  auto it = m_hitmap.find(key);
  if (it != m_hitmap.end())
  {
    return it->second;
  } else {
    return nullptr;
  }
}

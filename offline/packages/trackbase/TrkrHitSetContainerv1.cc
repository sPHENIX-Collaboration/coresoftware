/**
 * @file trackbase/TrkrHitSetContainerv1.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation for TrkrHitSetContainerv1
 */
#include "TrkrHitSetContainerv1.h"

#include "TrkrDefs.h"
#include "TrkrHitSetv1.h"

#include <cstdlib>

void TrkrHitSetContainerv1::Reset()
{
  for(auto&& [key, hitset] : m_hitmap)
    { delete hitset;}
  
  m_hitmap.clear();
}

void TrkrHitSetContainerv1::identify(std::ostream& os) const
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

TrkrHitSetContainerv1::ConstIterator
TrkrHitSetContainerv1::addHitSet(TrkrHitSet* newhit)
{ return addHitSetSpecifyKey(newhit->getHitSetKey(), newhit); }

TrkrHitSetContainerv1::ConstIterator
TrkrHitSetContainerv1::addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet* newhit)
{
  const auto ret = m_hitmap.insert(std::make_pair(key, newhit));
  if ( !ret.second )
  {
    std::cout << "TrkrHitSetContainerv1::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  } else {
    return ret.first;
  }
}

void TrkrHitSetContainerv1::removeHitSet(TrkrDefs::hitsetkey key)
{ 
  auto iter = m_hitmap.find( key );
  if(iter != m_hitmap.end()) 
    {
      TrkrHitSet* hitset = iter->second;
      delete hitset;
      m_hitmap.erase(iter); 
    }
}

void TrkrHitSetContainerv1::removeHitSet(TrkrHitSet *hitset)
{ removeHitSet( hitset->getHitSetKey() ); }

TrkrHitSetContainerv1::ConstRange
TrkrHitSetContainerv1::getHitSets(const TrkrDefs::TrkrId trackerid) const
{
  const TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid);
  const TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid);
  return std::make_pair( m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi) );
}

TrkrHitSetContainerv1::ConstRange
TrkrHitSetContainerv1::getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
{
  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid, layer);
  return std::make_pair( m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi) );
}

TrkrHitSetContainerv1::ConstRange
TrkrHitSetContainerv1::getHitSets() const
{ return std::make_pair(m_hitmap.cbegin(), m_hitmap.cend()); }

TrkrHitSetContainerv1::Iterator
TrkrHitSetContainerv1::findOrAddHitSet(TrkrDefs::hitsetkey key)
{
  auto it = m_hitmap.lower_bound( key );
  if( it == m_hitmap.end() || (key < it->first ) )
  {
    it = m_hitmap.insert(it, std::make_pair(key, new TrkrHitSetv1));
    it->second->setHitSetKey( key );
  }
  return it;
}

TrkrHitSet*
TrkrHitSetContainerv1::findHitSet(TrkrDefs::hitsetkey key)
{
  auto it = m_hitmap.find(key);
  if (it != m_hitmap.end())
  {
    return it->second;
  } else {
    return nullptr;
  }
}

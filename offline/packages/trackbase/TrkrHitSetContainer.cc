/**
 * @file trackbase/TrkrHitSetContainer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation for TrkrHitSetContainer
 */
#include "TrkrHitSetContainer.h"

#include "TrkrDefs.h"
#include "TrkrHitSet.h"

#include <cstdlib>

TrkrHitSetContainer::TrkrHitSetContainer()
{
}

TrkrHitSetContainer::~TrkrHitSetContainer()
{
  Reset();
}

void TrkrHitSetContainer::Reset()
{
  while (m_hitmap.begin() != m_hitmap.end())
  {
    delete m_hitmap.begin()->second;   // frees up memory for TrkrHit objects
    m_hitmap.erase(m_hitmap.begin());
  }
  return;
}

void TrkrHitSetContainer::identify(std::ostream& os) const
{
  ConstIterator iter;
  os << "Number of hits: " << size() << std::endl;
  for (iter = m_hitmap.begin(); iter != m_hitmap.end(); ++iter)
  {
    int layer = TrkrDefs::getLayer(iter->first);
    os << "hitsetkey " << iter->first << " layer " <<layer << std::endl;
    (iter->second)->identify();
  }
  return;
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSet(TrkrHitSet* newhit)
{
  return addHitSetSpecifyKey(newhit->getHitSetKey(), newhit);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet* newhit)
{
  auto ret = m_hitmap.insert(std::make_pair(key, newhit));
  if ( !ret.second )
  {
    std::cout << "TrkrHitSetContainer::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(const TrkrDefs::TrkrId trackerid) const
{
  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << std::endl;

  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid);

  ConstRange retpair;
  retpair.first = m_hitmap.lower_bound(keylo);
  retpair.second = m_hitmap.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(const TrkrDefs::TrkrId trackerid,
                                const char layer) const
{
  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = (tmp << TrackerDefs::bitshift_trackerid);
  // tmp = layer;
  // keylow |= (tmp << TrackerDefs::bitshift_layer);
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_layer) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << std::endl;

  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid, layer);

  ConstRange retpair;
  retpair.first = m_hitmap.lower_bound(keylo);
  retpair.second = m_hitmap.upper_bound(keyhi);
  return retpair;
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(void) const
{
  return std::make_pair(m_hitmap.begin(), m_hitmap.end());
}

TrkrHitSetContainer::Iterator
TrkrHitSetContainer::findOrAddHitSet(TrkrDefs::hitsetkey key)
{
  auto it = m_hitmap.lower_bound( key );
  if( it == m_hitmap.cend() || (key < it->first ) )
  {
    it = m_hitmap.insert(it, std::make_pair(key, new TrkrHitSet()));
    it->second->setHitSetKey( key );
  }
  return it;
}

TrkrHitSet*
TrkrHitSetContainer::findHitSet(TrkrDefs::hitsetkey key)
{
  TrkrHitSetContainer::ConstIterator it = m_hitmap.find(key);

  if (it != m_hitmap.end())
  {
    return it->second;
  }

  return 0;
}

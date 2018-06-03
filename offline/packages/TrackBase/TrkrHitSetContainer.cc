#include "TrkrHitSetContainer.h"

#include "TrkrDefUtil.h"
#include "TrkrHitSetv1.h"

#include <cstdlib>

TrkrHitSetContainer::TrkrHitSetContainer()
{
}

void TrkrHitSetContainer::Reset()
{
  while (m_hitmap.begin() != m_hitmap.end())
  {
    delete m_hitmap.begin()->second;
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
    os << "hit key 0x" << std::hex << iter->first << std::dec << std::endl;
    (iter->second)->identify();
  }
  return;
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSet(TrkrHitSet* newhit)
{
  TrkrDefs::hitsetkey key = newhit->getHitSetKey();
  if (m_hitmap.find(key) != m_hitmap.end())
  {
    TrkrDefUtil util;
    std::cout << "overwriting hit 0x" << std::hex << key << std::dec << std::endl;
    std::cout << "tracker id: " << util.getTrkrId(key) << std::endl;
  }
  m_hitmap[key] = newhit;
  return m_hitmap.find(key);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet* newhit)
{
  if (m_hitmap.find(key) != m_hitmap.end())
  {
    std::cout << "TrkrHitSetContainer::AddHitSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  newhit->setHitSetKey(key);
  m_hitmap[key] = newhit;
  return m_hitmap.find(key);
}

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(const TrkrDefs::TrkrId trackerid) const
{
  // TrkrDefs::hitsetkey tmp = trackerid;
  // TrkrDefs::hitsetkey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::hitsetkey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << std::endl;

  TrkrDefUtil util;
  TrkrDefs::hitsetkey keylo = util.getHitSetKeyLo(trackerid);
  TrkrDefs::hitsetkey keyhi = util.getHitSetKeyHi(trackerid);

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

  TrkrDefUtil util;
  TrkrDefs::hitsetkey keylo = util.getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = util.getHitSetKeyHi(trackerid, layer);

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
  TrkrHitSetContainer::Iterator it = m_hitmap.find(key);
  if (it == m_hitmap.end())
  {
    m_hitmap[key] = new TrkrHitSetv1();
    it = m_hitmap.find(key);
    TrkrHitSet* mhit = it->second;
    mhit->setHitSetKey(key);
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

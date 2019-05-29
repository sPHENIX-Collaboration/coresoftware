#ifndef TRACKBASE_TRKRHITSETCONTAINER_H
#define TRACKBASE_TRKRHITSETCONTAINER_H

#include "TrkrDefs.h"        // for hitsetkey, TrkrId

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

class TrkrHitSet;

/**
 * Container for TrkrHitSet objects
 */
class TrkrHitSetContainer : public PHObject
{
 public:
  typedef std::map<TrkrDefs::hitsetkey, TrkrHitSet *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  //! ctor
  TrkrHitSetContainer();

  //! dtor
  virtual ~TrkrHitSetContainer();
  //! PHObject functions
  void Reset();
  void identify(std::ostream &os = std::cout) const;

  //! Add a TrkrHitSet to the container
  ConstIterator addHitSet(TrkrHitSet *newHit);
  ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet *newHit);

  //! preferred removal method, key is currently the hit id
  void removeHitSet(TrkrDefs::hitsetkey key)
  {
    m_hitmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void removeHitSet(TrkrHitSet *hit)
  {
    Iterator its = m_hitmap.begin();
    Iterator last = m_hitmap.end();
    for (; its != last;)
    {
      if (its->second == hit)
      {
        m_hitmap.erase(its++);
      }
      else
      {
        ++its;
      }
    }
  }

  //! find or add HitSet
  Iterator findOrAddHitSet(TrkrDefs::hitsetkey key);

  //! return all HitSets matching a given detid
  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const;

  //! return all HitSets matching a given detid, layer
  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const char layer) const;

  //! return all HitSets
  ConstRange getHitSets(void) const;

  //! return a given HitSet based on its key
  TrkrHitSet *findHitSet(TrkrDefs::hitsetkey key);

  unsigned int size(void) const
  {
    return m_hitmap.size();
  }

 protected:
  Map m_hitmap;
  ClassDef(TrkrHitSetContainer, 1)
};

#endif //TRACKBASE_TRKRHITSETCONTAINER_H

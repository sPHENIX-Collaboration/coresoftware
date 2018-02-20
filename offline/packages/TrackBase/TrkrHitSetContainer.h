#ifndef __TrkrHitSetContainer_H__
#define __TrkrHitSetContainer_H__

#include "TrkrHitSet.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

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
  virtual ~TrkrHitSetContainer() {}
  //! PHObject functions
  void Reset();
  void identify(std::ostream &os = std::cout) const;

  //! Add a TrkrHitSet to the container
  ConstIterator AddHitSet(TrkrHitSet *newHit);
  ConstIterator AddHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet *newHit);

  //! preferred removal method, key is currently the hit id
  void RemoveHitSet(TrkrDefs::hitsetkey key)
  {
    hitmap_.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveHitSet(TrkrHitSet *hit)
  {
    Iterator its = hitmap_.begin();
    Iterator last = hitmap_.end();
    for (; its != last;)
    {
      if (its->second == hit)
      {
        hitmap_.erase(its++);
      }
      else
      {
        ++its;
      }
    }
  }

  //! find or add HitSet
  Iterator FindOrAddHitSet(TrkrDefs::hitsetkey key);

  //! return all HitSets matching a given detid
  ConstRange GetHitSets(const TrkrDefs::TRKRID trackerid) const;

  //! return all HitSets matching a given detid, layer
  ConstRange GetHitSets(const TrkrDefs::TRKRID trackerid, const char layer) const;

  //! return all HitSets
  ConstRange GetHitSets(void) const;

  //! return a given HitSet based on its key
  TrkrHitSet *FindHitSet(TrkrDefs::hitsetkey key);

  unsigned int size(void) const
  {
    return hitmap_.size();
  }

 protected:
  Map hitmap_;
  ClassDef(TrkrHitSetContainer, 1)
};

#endif

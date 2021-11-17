#ifndef TRACKBASE_TRKRHITSETCONTAINER_H
#define TRACKBASE_TRKRHITSETCONTAINER_H
/**
 * @file trackbase/TrkrHitSetContainerv1.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * base class for hitset container
 */

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
  
  using Map = std::map<TrkrDefs::hitsetkey, TrkrHitSet *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  //! dtir
  ~TrkrHitSetContainer() override = default;

  //! PHObject functions
  void Reset()  override;
  
  //! Add a TrkrHitSet to the container
  virtual ConstIterator addHitSet(TrkrHitSet*);

  virtual ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, TrkrHitSet*);

  //! preferred removal method, key is currently the hit id
  virtual void removeHitSet(TrkrDefs::hitsetkey)
  {}

  //! inefficent, use key where possible instead
  virtual void removeHitSet(TrkrHitSet*)
  {}

  //! find or add HitSet
  virtual Iterator findOrAddHitSet(TrkrDefs::hitsetkey);

  //! return all HitSets matching a given detid
  virtual ConstRange getHitSets(const TrkrDefs::TrkrId) const;

  //! return all HitSets matching a given detid, layer
  virtual ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const;

  //! return all HitSets
  virtual ConstRange getHitSets() const;

  //! return a given HitSet based on its key
  virtual TrkrHitSet *findHitSet(TrkrDefs::hitsetkey)
  { return nullptr; }

  virtual unsigned int size() const
  { return 0; }

  protected:
  //! ctor
  TrkrHitSetContainer() = default;

  private:

  ClassDefOverride(TrkrHitSetContainer, 1)
};

#endif //TRACKBASE_TRKRHITSETCONTAINER_H

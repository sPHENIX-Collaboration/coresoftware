#ifndef TRACKBASE_RAWHITSETCONTAINER_H
#define TRACKBASE_RAWHITSETCONTAINER_H
/**
 * @file trackbase/RawHitSetContainerv1.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * base class for hitset container
 */

#include "TrkrDefs.h"        // for hitsetkey, TrkrId

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

class RawHitSet;

/**
 * Container for RawHitSet objects
 */
class RawHitSetContainer : public PHObject
{
 public:
  
  using Map = std::map<TrkrDefs::hitsetkey, RawHitSet *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  //! dtir
  ~RawHitSetContainer() override = default;

  //! PHObject functions
  void Reset()  override;
  
  //! Add a RawHitSet to the container
  virtual ConstIterator addHitSet(RawHitSet*);

  virtual ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, RawHitSet*);

  //! preferred removal method, key is currently the hit id
  virtual void removeHitSet(TrkrDefs::hitsetkey)
  {}

  //! inefficent, use key where possible instead
  virtual void removeHitSet(RawHitSet*)
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
  virtual RawHitSet *findHitSet(TrkrDefs::hitsetkey)
  { return nullptr; }

  virtual unsigned int size() const
  { return 0; }

  protected:
  //! ctor
  RawHitSetContainer() = default;

  private:

  ClassDefOverride(RawHitSetContainer, 1)
};

#endif //TRACKBASE_RAWHITSETCONTAINER_H

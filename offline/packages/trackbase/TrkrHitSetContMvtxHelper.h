/**
 * @file trackbase/TrkrHitSetContMvtxHelper.h
 * @author Yasser Corrales Morales <ycmorales@bnl.gov>
 * @date Febraury 2025
 * base class for Mvtx hitsetkey container per strobe
 */

#ifndef TRACKBASE_TRKRHITSETCONTMVTXHELPER_H
#define TRACKBASE_TRKRHITSETCONTMVTXHELPER_H

#include "TrkrDefs.h"  // for hitsetkey, TrkrId

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>       // for map
#include <memory>
#include <set>      // for set
#include <utility>  // for pair

// class TrkrHitSetKey;

/**
 * Container for Mvtx hitsekey per strobe
 */
class TrkrHitSetContMvtxHelper : public PHObject
{
 public:
  using map_tp_value = std::set<TrkrDefs::hitsetkey>;
  using Map = std::map<int32_t, map_tp_value>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;

  //! cp/mv ctor
  TrkrHitSetContMvtxHelper(const TrkrHitSetContMvtxHelper &) = default;
  TrkrHitSetContMvtxHelper(TrkrHitSetContMvtxHelper &&) = default;
  //! cp/mv assignment
  TrkrHitSetContMvtxHelper &operator=(const TrkrHitSetContMvtxHelper &) = default;
  TrkrHitSetContMvtxHelper &operator=(TrkrHitSetContMvtxHelper &&) = default;

  //! dtor
  ~TrkrHitSetContMvtxHelper() override = default;

  //! PHObject functions
  void Reset() override;

  virtual bool addHitSetKey(const int32_t & /*dummy*/, const TrkrDefs::hitsetkey & /*dummy*/)
  {
    return false;
  }

  //! preferred removal method, key is currently the hit id
  virtual bool removeHitSetKey(int32_t /*strobe*/, const TrkrDefs::hitsetkey /*dummy*/)
  {
    return false;
  }

  //! return all HitSetKeys
  virtual const map_tp_value &getHitSetKeys(int32_t /*strobe*/) = 0;

  virtual unsigned int size() const
  {
    return 0;
  }

 protected:
  //! ctor
  TrkrHitSetContMvtxHelper() = default;

 private:
  ClassDefOverride(TrkrHitSetContMvtxHelper, 1)
};

#endif  // TRACKBASE_TRKRHITSETCONTAINER_H

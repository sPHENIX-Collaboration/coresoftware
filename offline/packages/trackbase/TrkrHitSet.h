#ifndef TRACKBASE_TRKRHITSET_H
#define TRACKBASE_TRKRHITSET_H
/**
 * @file trackbase/TrkrHitSet.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date 4 June 2018
 * @brief Base Class Container for storing TrkrHit's
 */

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>  // for pair

//! forward declaration
class TrkrHit;

/**
 * @brief Container for storing TrkrHit's
 *
 * Container object which stores a set of TrkrHit objects.
 * A single TrkrHitSet is meant to represent a geometric detector object
 * which bounds clustering. Therefore, a TrkrHitSet should contain all
 * TrkrHits which could belong the the same cluster.
 */
class TrkrHitSet : public PHObject
{
 public:
  // iterator typedef
  using Map = std::map<TrkrDefs::hitkey, TrkrHit*>;
  using ConstIterator = Map::const_iterator;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  //! TObject functions
  void identify(std::ostream& /*os*/ = std::cout) const override
  {
  }

  void Reset() override
  {
  }

  /**
   * @brief Set the key for this object
   * @param key
   */
  virtual void setHitSetKey(const TrkrDefs::hitsetkey)
  {
  }

  /**
   * @brief Get the key for this object
   * @param[out] object key
   */
  virtual TrkrDefs::hitsetkey getHitSetKey() const
  {
    return TrkrDefs::HITSETKEYMAX;
  }

  /**
   * @brief Add a hit to this container using a specific key.
   * @param[in] key Hit key
   * @param[in] hit Hit to be added.
   *
   * NOTE: This TrkrHitSet takes ownership of the passed TrkrHit pointer
   * and will delete it in the Reset() method.
   */
  virtual ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*);

  /**
   * @brief Remove a hit using its key
   * @param[in] key to be removed
   */
  virtual void removeHit(TrkrDefs::hitkey)
  {
  }

  /**
   * @brief Get a specific hit based on its index.
   * @param key of the desired hit
   * @param[out] Pointer to the desired hit. nullptr if no hit.
   *
   * Get a desired hit based on its key.
   */
  virtual TrkrHit* getHit(const TrkrDefs::hitkey) const
  {
    return nullptr;
  }

  /**
   * @brief Get all hits
   * @param[out] Pair of iterator to vector begin and end
   */
  virtual ConstRange getHits() const;

  /**
   * @brief Get the number of hits stored
   * @param[out] number of hits
   */
  // Get number of hits
  virtual unsigned int size() const
  {
    return 0;
  }

  protected:

  //! ctor, not to be called
  TrkrHitSet() = default;

private:
  ClassDefOverride(TrkrHitSet, 1);
};

#endif  //TRACKBASE_TRKRHITSET_H

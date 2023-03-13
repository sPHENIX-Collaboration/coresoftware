#ifndef TRACKBASE_RAWHITSET_H
#define TRACKBASE_RAWHITSET_H
/**
 * @file trackbase/RawHitSet.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date 4 June 2018
 * @brief Base Class Container for storing RawHit's
 */

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>  // for pair

//! forward declaration
class RawHit;
class RawHitTpc;

/**
 * @brief Container for storing RawHit's
 *
 * Container object which stores a set of RawHit objects.
 * A single RawHitSet is meant to represent a geometric detector object
 * which bounds clustering. Therefore, a RawHitSet should contain all
 * RawHits which could belong the the same cluster.
 */
class RawHitSet : public PHObject
{
 public:
  // iterator typedef
  using Vector = std::vector<RawHit*>;
  using VectorTpc2D = std::vector<std::vector<uint8_t>>;
  using ConstIterator = Vector::const_iterator;
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
   * NOTE: This RawHitSet takes ownership of the passed RawHit pointer
   * and will delete it in the Reset() method.
   */
  virtual void addHit(RawHit*);
  //  virtual void addTpcHit(unsigned short phibin,RawHit*);
  virtual void setTpcPhiBins(unsigned short phibins);
  /**
   * @brief Remove a hit using its key
   * @param[in] key to be removed
   */
  /* virtual void removeHit(TrkrDefs::hitkey)
  {
  }
  */
  /**
   * @brief Get a specific hit based on its index.
   * @param key of the desired hit
   * @param[out] Pointer to the desired hit. nullptr if no hit.
   *
   * Get a desired hit based on its key.
   */
  /*  virtual RawHit* getHit(const TrkrDefs::hitkey) const
  {
    return nullptr;
  }
  */
  /**
   * @brief Get all hits
   * @param[out] Pair of iterator to vector begin and end
   */
  virtual ConstRange getHits() const;
  //  virtual ConstRange getTpcHits(unsigned short phibin) const;
  /**
   * @brief Get the number of hits stored
   * @param[out] number of hits
   */
  // Get number of hits
  virtual unsigned int size() const
  {
    return 0;
  }
  virtual unsigned int tpcphibins() const
  {
    return 0;
  }
  protected:

  //! ctor, not to be called
  RawHitSet() = default;

private:
  ClassDefOverride(RawHitSet, 1);
};

#endif  //TRACKBASE_RAWHITSET_H

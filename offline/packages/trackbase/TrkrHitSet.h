/**
 * @file trackbase/TrkrHitSet.h
 * @author D. McGlinchey
 * @date 4 June 2018
 * @brief Container for storing TrkrHit's
 */
#ifndef TRACKBASE_TRKRHITSET_H
#define TRACKBASE_TRKRHITSET_H

#include "TrkrDefUtil.h"

#include <TObject.h>

#include <iostream>
#include <vector>

class TrkrHit;

/**
 * @brief Container for storing TrkrHit's
 *
 * Container object which stores a set of TrkrHit objects.
 * A single TrkrHitSet is meant to represent a geometric detector object
 * which bounds clustering. Therefore, a TrkrHitSet should contain all
 * TrkrHits which could belong the the same cluster.
 */
class TrkrHitSet : public TObject
{
 public:
  // iterator typedef
  typedef std::vector<TrkrHit*> Vec;
  typedef Vec::const_iterator ConstIterator;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  
  //! ctor
  TrkrHitSet();
  //! dtor
  virtual ~TrkrHitSet() {}
  //! TObject functions
  void identify(std::ostream& os = std::cout) const;
  void Reset();
  void print() const;

  /**
   * @brief Set the key for this object
   * @param key
   */
  void setHitSetKey(const TrkrDefs::hitsetkey key) { m_hitSetKey = key; }

  /**
   * @brief Get the key for this object
   * @param[out] object key
   */
  TrkrDefs::hitsetkey getHitSetKey() const { return m_hitSetKey; }

  /**
   * @brief Add a hit to this container.
   * @param hit to be added.
   * @param[out] index of the added hit, which can be used for access.
   *
   * NOTE: This TrkrHitSet takes ownership of the passed TrkrHit pointer
   * and will delete it in the Reset() method.
   */
  unsigned int addHit(TrkrHit* hit);

  /**
   * @brief Get a specific hit based on its index.
   * @param index of the desired hit
   * @param[out] Pointer to the desired hit.
   *
   * Get a desired hit based on its index in the internal 
   * storage vector.
   */
  TrkrHit* getHit(unsigned int ihit);
  
  /**
   * @brief Get all hits
   * @param[out] Pair of iterator to vector begin and end
   */
  ConstRange getHits();

  /**
   * @brief Get the number of hits stored
   * @param[out] number of hits
   */
  // Get number of hits
  unsigned int size() { return m_hits.size(); }
   
 private:
  TrkrDefs::hitsetkey m_hitSetKey; /// unique key for this object
  Vec m_hits; /// storage for TrkrHit objects
  ClassDef(TrkrHitSet, 1);
};

#endif //TRACKBASE_TRKRHITSET_H

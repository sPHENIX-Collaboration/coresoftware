/**
 * @file trackbase/TrkrClusterHitAssoc.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Base class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERHITASSOC_H
#define TRACKBASE_TRKRCLUSTERHITASSOC_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair
#include <climits>

/**
 * @brief Base class for associating clusters to the hits that went into them
 *
 * Store the associations between clusters and the hits that went into them.
 */
class TrkrClusterHitAssoc : public PHObject
{
public:

  using Map = std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>;
  using ConstIterator = Map::const_iterator;
  using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;
  
  TrkrClusterHitAssoc() = default;

  virtual void Reset() = 0;

  virtual void identify(std::ostream &os = std::cout) const = 0;
  
  /**
   * @brief Add association between cluster and hit
   * @param[in] ckey Cluster key
   * @param[in] hidx Index of the hit in TrkrHitSet
   */
  virtual void addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx) = 0;

  //! get pointer to cluster-to-hit map corresponding to a given hitset id
  virtual Map* getClusterMap(TrkrDefs::hitsetkey) {return nullptr;}

  /**
   * @brief Get all the hits associated with a cluster by key
   * @param[in] ckey Cluster key
   * @param[out] Range over hits associated with @c ckey
   */

  virtual ConstRange getHits(TrkrDefs::cluskey) = 0;

  virtual unsigned int size() const {return 0;}

private:

  ClassDef(TrkrClusterHitAssoc, 1);
};

#endif // TRACKBASE_TRKRCLUSTERHITASSOC_H

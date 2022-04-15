/**
 * @file trackbase/TrkrClusterCrossingAssoc.h
 * @author Tony Frawley
 * @date March 2022
 * @brief Base class for associating clusters to the bunch crossing they originated from
 */
#ifndef TRACKBASE_TRKRCLUSTERCROSSINGASSOC_H
#define TRACKBASE_TRKRCLUSTERCROSSINGASSOC_H

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
class TrkrClusterCrossingAssoc : public PHObject
{
public:

  using Map = std::multimap<TrkrDefs::cluskey, short int>;
  using ConstIterator = Map::const_iterator;
  using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;
  
  void Reset() override;

  /**
   * @brief Add association between cluster and hit
   * @param[in] ckey Cluster key
   * @param[in] bunch crossing number
   */
  virtual void addAssoc(TrkrDefs::cluskey ckey, short int hidx) = 0;

  virtual ConstRange getCrossings(TrkrDefs::cluskey) const = 0;

  virtual ConstRange getAll() const = 0;

  virtual unsigned int size() const {return 0;}

protected:
  TrkrClusterCrossingAssoc() = default;

private:

  ClassDefOverride(TrkrClusterCrossingAssoc, 1);
};

#endif // TRACKBASE_TRKRCLUSTERCROSSINGASSOC_H

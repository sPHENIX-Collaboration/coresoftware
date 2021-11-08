/**
 * @file trackbase/TrkrClusterIterationMap.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Base class for associating clusters to the hits that went into them
 */
#ifndef TRACKBASE_TRKRCLUSTERITERATIONMAP_H
#define TRACKBASE_TRKRCLUSTERITERATIONMAP_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair
#include <climits>

/**
 * @brief Base class for associating clusters to iterations they were used in
 *
 * Store the associations between clusters and trackign iterations.
 */
class TrkrClusterIterationMap : public PHObject
{
public:

  using Map = std::map<TrkrDefs::cluskey, short int>;
  //  using ConstIterator = Map::const_iterator;
  // using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;
  
  void Reset() override;

  /**
   * @brief Add association between cluster andthe tracking iteration it was used in
   * @param[in] ckey Cluster key
   * @param[in] tracking iteration
   */
  virtual void addIteration(TrkrDefs::cluskey ckey, short int iter) = 0;

  virtual short int getIteration(TrkrDefs::cluskey ckey)  = 0;

  virtual unsigned int size() const {return 0;}

protected:
  TrkrClusterIterationMap() = default;

private:

  ClassDefOverride(TrkrClusterIterationMap, 1);
};

#endif // TRACKBASE_TRKRCLUSTERITERATIONMAP_H

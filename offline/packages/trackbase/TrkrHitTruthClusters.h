#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERS_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERS_H

/**
 * @file trackbase/TrkrHitCellAssoc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Association object for PHG4Hits contributiong to TrkrHits
 */

#include "TrkrDefs.h"

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Clusters contributiong to TrkrHitTruthClusters
 *
 * Association object holding a vector of PHG4Clusters associated with given Truth Electrons
 * and associated streamed electrons
 *
 */
class TrkrHitTruthClusters : public PHObject
{
  
  public:
  
  //! typedefs for convenience 
  using MMap = std::map< int /*track-id*/, 
                         std::vector<std::array<float,6>> /*array of mean-std of phi,eta,z*/>;
  using Iterator = MMap::iterator;
  using ConstIterator = MMap::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  void Reset() override
  {}

  virtual void print_clusters (std::ostream &/*os*/ = std::cout) const 
  {}

  virtual void push_truth_cluster(const int /*track_id*/, 
          const std::array<double,6>& /*phi_eat_z_data*/, const double /*sum_E*/) 
  {}

  //virtual void removeAssoc(const TrkrDefs::hitsetkey /*hitsetkey*/, const TrkrDefs::hitkey /*hitkey*/)
  //{}

  protected:
  //! ctor
  TrkrHitTruthClusters() = default;

  private:

  ClassDefOverride(TrkrHitTruthClusters, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H

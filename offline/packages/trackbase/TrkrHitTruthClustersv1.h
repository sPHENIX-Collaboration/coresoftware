#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H

/**
 * @file trackbase/TrkrHitTruthClustersv1.h
 * @author D. STEWART
 * @date April 2022
 */

#include "TrkrDefs.h"
#include "TrkrHitTruthClusters.h"

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

/**
 * @brief Association object for PHG4Cells contributiong to TrkrHits
 *
 * Association object holding a multimap of PHG4Cells associated with a given TrkrHit
 */
class TrkrHitTruthClustersv1 : public TrkrHitTruthClusters
{
  
  public:

  TrkrHitTruthClustersv1() = default;
 
  void Reset() override;

  void print_clusters (std::ostream &os = std::cout) const override;

  /* CentroidsFor1Track& add_track_centroids(const int track_id) override; */
  /* void add_track_centroids(short track_id, VecEC& centroids) override; */
  VecEC& get_new_centroids_vec (short track_id) override;

  MMap  m_map;

  private:
  ClassDefOverride(TrkrHitTruthClustersv1, 1);

};

#endif //TRACKBASE_TRKRHITTRUTHCLUSTERS_H


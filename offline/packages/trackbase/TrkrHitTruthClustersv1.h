#ifndef TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H
#define TRACKBASE_TRKRHITTRUTHCLUSTERSV1_H

/**
 * @file trackbase/TrkrHitTruthClustersv1.h
 * @author D. STEWART
 * @date April 2022
 */

#include "TrkrHitTruthClusters.h"

#include "TrkrDefs.h"

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair

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
  void identify(std::ostream& os = std::cout) const override {print_clusters(os);}

  void print_clusters(std::ostream& os = std::cout) const override;

  void push_truth_cluster(const int track_id,
                          const std::array<double, 8>& phi_eta_z_data, const double sum_E) override;

 private:
  MMap m_map;

  ClassDefOverride(TrkrHitTruthClustersv1, 1);
};

#endif  //TRACKBASE_TRKRHITTRUTHCLUSTERS_H

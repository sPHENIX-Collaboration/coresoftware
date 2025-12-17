#ifndef TRACKCLUSEVALUATOR_H
#define TRACKCLUSEVALUATOR_H
// A class that will collect the clusters for SVTX and PHG$ tracks and compare them.
// The actual comparison, cluster to cluster, is from a pointer to a user provided TrkrClusterIsMatcher
#include "TrkrClusLoc.h"

#include <trackbase/TrkrDefs.h>

#include <array>
#include <vector>

class TrkrClusterIsMatcher;
class SvtxTrack;
class TrkrTruthTrack;
class TrkrClusterContainer;

class TrackClusEvaluator
{
 private:
  using Vector = std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::cluskey>>;

  using Iter = Vector::iterator;

  /* TrkrClusterComparer*  comp      {nullptr}; // DEPRECATED */

  std::array<int, 5> cntclus(Vector& keys);                                      // nclusters MVTX, INTT, TPC, TPOT, TOTAL
  std::array<int, 5> cnt_matchedclus(Vector& keys, std::vector<bool>& matches);  // same but number matched

 public:
  TrkrClusterIsMatcher* ismatcher{nullptr};
  /* void set_ismatcher(TrkrClusterIsMatcher* _ismatcher) { ismatcher = _ismatcher; }; */
  TrackClusEvaluator(TrkrClusterIsMatcher* tcm = nullptr)
    : ismatcher{tcm} {};

  TrkrClusterContainer* get_PHG4_clusters() const;
  TrkrClusterContainer* get_SVTX_clusters() const;

  Vector svtx_keys{};
  Vector phg4_keys{};

  bool collect_match_statistic = false;
  double match_stat{0};

  void reset();
  std::array<int, 3> find_matches();  // populated matches_{svtx,phg4};
                                      // return's {n-matched, n-phg4, n-svtx}
  std::array<int, 3> find_matches(TrkrTruthTrack* g4_track, SvtxTrack* sv_track);

  int phg4_n_matched();  // also same as phg4_cnt_matchedclus()[4]
  int svtx_n_matched();  // should be almost always the same
                         // which is ALMOST guaranteed to be same as svtx_cnt_matchedclus()[4]
  int phg4_nclus() { return (int) phg4_keys.size(); }
  int svtx_nclus() { return (int) svtx_keys.size(); }

  std::vector<bool> svtx_matches;
  std::vector<bool> phg4_matches;

  int addClusKeys(SvtxTrack*);       // return number of clusters
  int addClusKeys(TrkrTruthTrack*);  // return number of clusters

  std::array<int, 5> svtx_cntclus() { return cntclus(svtx_keys); };  // Mvtx Intt Tpc TPOT Sum
  std::array<int, 5> phg4_cntclus() { return cntclus(phg4_keys); };

  std::array<int, 5> svtx_cnt_matchedclus() { return cnt_matchedclus(svtx_keys, svtx_matches); };
  std::array<int, 5> phg4_cnt_matchedclus() { return cnt_matchedclus(phg4_keys, phg4_matches); };

  // Convenience functions, asking for cluster locations
  std::vector<TrkrClusLoc> phg4_clusloc_all();
  std::vector<TrkrClusLoc> phg4_clusloc_unmatched();
  std::vector<TrkrClusLoc> svtx_clusloc_all();
  std::vector<TrkrClusLoc> svtx_clusloc_unmatched();
  std::vector<TrkrClusLoc> clusloc_matched();
};

#endif

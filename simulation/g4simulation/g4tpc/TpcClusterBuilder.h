#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

// Revised 04-Dec-2022, David Stewart
// basic use:
//  (a) optionally used set_current_track() (pointless to not set the track though,
//      otherwise the new TrkrClsuters won't be assigned to the proper location)
//  (b) fill with TrkrHits in addhitset() (in the MapToPadPlane module)
//  (c) after all TrkrHit's have been added for a given truth track,
//      call cluster_and_reset(false) to generate the TrkrClusters from the
//      TrkrHit's, fill in the clusterkeys in the truth_track, and reset
//      the TrkrHit container.
//  (d) at the end of the event, run cluster_and_reset(true) which will check
//      if there are clusters to make, and also clear out the counter for the
//      hitsetkeys
//
//  Note:
//  - the algorithm to get the local surface will probably be updated in
//    TpcClusterizer, when that happens update here as well

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainerv1.h>

#include <phool/PHObject.h>

#include <climits>
#include <iostream>
#include <map>

class ClusHitsVerbosev1;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrTruthTrack;
class TrkrTruthTrackContainer;

// This is the basic data for each set of TrkrHits from each TrkrHitsSet
// to be used in tpc/TpcClusterizer.cc
class TpcClusterBuilder
{
 public:
  TpcClusterBuilder(){};
  ~TpcClusterBuilder()
  {
    delete m_hits;
  };

  void fixme_check();
  void fixme_short_check();

  bool b_collect_hits{false};
  bool needs_input_nodes() { return m_needs_input_nodes; }
  void cluster_hits(TrkrTruthTrack* track);
  void addhitset(TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void set_current_track(TrkrTruthTrack* _trkrtruthtrack);
  void print(TrkrTruthTrackContainer*, int nclusprint = -1);
  void print_file(TrkrTruthTrackContainer*, const std::string&);
  void set_verbosity(int verbosity_level) { verbosity = verbosity_level; }

  void clear_hitsetkey_cnt();
  void set_pixel_thresholdrat(double val) { m_pixel_thresholdrat = val; };
  void set_input_nodes(TrkrClusterContainer* _truth_cluster_container, ActsGeometry* _ActsGeometry, PHG4TpcCylinderGeomContainer* _geom_container, ClusHitsVerbosev1* _clushitsverbose);

 private:
  ActsGeometry* m_tGeometry{nullptr};  // used to generate clusters
  ClusHitsVerbosev1* mClusHitsVerbose{nullptr};
  PHG4TpcCylinderGeomContainer* geom_container{nullptr};
  TrkrClusterContainer* m_clusterlist{nullptr};  // fill for output

  // internal containers to fill and consume hits and fill with tracks
  TrkrHitSetContainer* m_hits{new TrkrHitSetContainerv1()};
  /* TrkrTruthTrack*      current_track { nullptr }; */
  std::map<TrkrDefs::hitsetkey, unsigned int> hitsetkey_cnt{};

  int n_tracks{0};
  int verbosity{0};

  double AdcClockPeriod{53.0};  // ns

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  double m_sampa_tbias{39.6};  // ns

  // for pixel thresholds
  double m_pixel_thresholdrat{0.01};

  bool m_needs_input_nodes{true};
};

#endif  // TRACKBASE_PADPLANEREADOUTSTRUCT_H

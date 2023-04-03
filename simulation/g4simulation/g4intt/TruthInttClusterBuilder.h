#ifndef G4INT_TRUTHINTTCLUSTERBUILDER
#define G4INT_TRUTHINTTCLUSTERBUILDER

// Generated March 2023, David Stewart
//
//     This file is to PHG4InttHitReco.cc & InttClusterizer.cc 
//     as g4tpc/TpcClusterBuilder.cc is to PHG4TpcElectronDrift.cc & TpcClusterizer.cc
//     and g4intt/TruthMvtxClusterBuilder.cc is to PHG4MvtxHitReco.cc & MvtxClusterizer.cc
//
//     Namely, it is a builder that can collect the TrkrHits generated in PHG4MvtxHitReco.cc
//     for each truth track, and then use the logic of MvtxClusterizer.cc to cluster them 
//     into TrkrTruthTracks and fill them in the TrkrTruthTrackContainer on the node tree.
//
// Purpose and use:
//    - Collect only TrkrHits in PHG4InttHitReco associated with truth tracks (one track at a time)
//    - Everytime there is a new truth track (or event track is done) process the existing hits from
//      the last truth track into new TrkrClusters associated with a newly generated TrkrTruthTrack
//      which is, in tern, added to a TrkrTruthTrackContainer.
//
//      The logic to generate the TrkrClusters is meant to be identical to that in offline/packages/intt/InttClusterizer.cc,
//      which generates the clusters for *all* TrkrHits, not just those associated with truth tracks.
//
//      This is similar to how the g4simulation/g4tpc/TpcClusterBuilder module works with PHG4TpcElectronDrift.cc
//      and packages/tpc/TpcClusterizer.cc

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <climits>
#include <map>
#include <vector>
#include <string>
#include <utility>

// input
class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class PHG4CylinderGeomContainer;
// !! class TrkrClusterHitAssoc;
// !! class TrkrClusterCrossingAssoc;
class TrkrHit;
/* class RawHit; */
/* class RawHitSet; */
/* class RawHitSetContainer; */

// output
class TrkrTruthTrackContainer;
class TrkrTruthTrack;
class PHG4TruthInfoContainer;
class PHG4Hit;

class TruthInttClusterBuilder {

  // member data
  int                      trkid {-1};
  bool                     is_emb {false};
  PHG4TruthInfoContainer  *m_truthinfo {nullptr};

  // intermediary hits, used to cluster Truth hits
  TrkrHitSetContainer     *m_hits { nullptr }; // save the hits as they come along
  std::map<TrkrDefs::hitsetkey,unsigned int> hitsetkey_cnt {}; // counter for filling the truthtracks

  //output nodes to fill
  TrkrClusterContainer    *m_clusters;
  TrkrTruthTrackContainer *m_truthtracks;
  /* PHG4TpcCylinderGeomContainer* m_geom_container */
  PHG4CylinderGeomContainer* m_geom_container { nullptr };

  // member data required for clustering

  public:
  TruthInttClusterBuilder(TrkrClusterContainer* _clusters, 
      TrkrTruthTrackContainer* _truth_tracks,
      int _verbosity );
  ~TruthInttClusterBuilder();

  void set_truthinfo (PHG4TruthInfoContainer *_truthinfo) { m_truthinfo = _truthinfo; };
  void set_geom_container (PHG4CylinderGeomContainer *_geom_container) { m_geom_container = _geom_container; };
  void check_g4hit   (PHG4Hit                *hit);
  void addhitset     (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void reset();

  void set_verbosity(int _verbosity) { m_verbosity = _verbosity; };

  private:
  void cluster_hits();
  int m_verbosity {0};
  

  bool ladder_are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer);
  /* bool ladder_are_adjacent(RawHit* lhs,  RawHit* rhs, const int layer); */

  void set_z_clustering(const int layer, const bool make_z_clustering)
  {
    _make_z_clustering.insert(std::make_pair(layer, make_z_clustering));
  }
  bool get_z_clustering(const int layer) const
  {
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) return true;
    return _make_z_clustering.find(layer)->second;
  }

  // settings
  /* float _fraction_of_mip; */
  std::map<int, float> _thresholds_by_layer;  // layer->threshold
  std::map<int, bool> _make_z_clustering;     // layer->z_clustering_option
  std::map<int, bool> _make_e_weights;        // layer->energy_weighting_option
  /* bool do_hit_assoc = true; */
  /* bool do_read_raw = false; */
  /* int m_cluster_version = 4; */
};

#endif

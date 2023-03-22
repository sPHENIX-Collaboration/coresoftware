#ifndef G4MVTX_TRUTHMVTXCLUSTERBUILDER__H
#define G4MVTX_TRUTHMVTXCLUSTERBUILDER__H

// This file is to PHG4MvtxHitReco.cc & MvtxClusterizer.cc 
// as g4tpc/TpcClusterBuilder.cc is to PHG4TpcElectronDrift.cc & TpcClusterizer.cc
// and g4intt/TruthInttClusterBuilder.cc is to PHG4InttHitReco.cc & InttClusterizer.cc
//
// Namely, it is a builder that can collect the TrkrHits generated in PHG4MvtxHitReco.cc
// for each truth track, and then use the logic of MvtxClusterizer.cc to cluster them 
// into TrkrTruthTracks and fill them in the TrkrTruthTrackContainer on the node tree.
//
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
class PHG4CylinderGeomContainer;
class TrkrClusterContainer;
class TrkrHitSetContainer;
// !! class TrkrClusterHitAssoc;
// !! class TrkrClusterCrossingAssoc;
class TrkrHit;
class RawHit;
class RawHitSet;
class RawHitSetContainer;

// output
class TrkrTruthTrackContainer;
class TrkrTruthTrack;
class PHG4TruthInfoContainer;
class PHG4Hit;

class MvtxHitPruner;

class TruthMvtxClusterBuilder {

  // member data
  int                      trkid {-1};
  bool                     is_emb {false};
  PHG4TruthInfoContainer  *m_truthinfo {nullptr};
  MvtxHitPruner *m_MvtxHitPruner { nullptr };

  // intermediary hits, used to cluster Truth hits
  TrkrHitSetContainer     *m_hits {}; // save the hits as they come along
  std::map<TrkrDefs::hitsetkey,unsigned int> hitsetkey_cnt {}; // counter for filling the truthtracks

  //output nodes to fill
  TrkrClusterContainer    *m_clusters;
  TrkrTruthTrackContainer *m_truthtracks;

  // member data required for clustering

  public:
  TruthMvtxClusterBuilder(TrkrClusterContainer* _clusters, 
      TrkrTruthTrackContainer* _truth_tracks,
      int _verbosity);
  ~TruthMvtxClusterBuilder(); 

  void set_truthinfo (PHG4TruthInfoContainer *_truthinfo) { m_truthinfo = _truthinfo; };
  void check_g4hit   (PHG4Hit                *hit);
  void addhitset     (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void reset();

  void set_verbosity(int _verbosity) { m_verbosity = _verbosity; };
  void set_geom_container (PHG4CylinderGeomContainer *_geom_container) { m_geom_container = _geom_container; };
  void set_HitPruner(MvtxHitPruner* _) { m_MvtxHitPruner = _; };

  void print_clusters(int nclusters=-1);

  //! option to turn off z-dimension clustering
  void SetZClustering(const bool make_z_clustering)
  {
    m_makeZClustering = make_z_clustering;
  }
  bool GetZClustering() const
  {
    return m_makeZClustering;
  }

  private:
  void cluster_hits();
  int m_verbosity {0};
  /* bool are_adjacent(RawHit* lhs, RawHit* rhs); */
  bool are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs);
  PHG4CylinderGeomContainer* m_geom_container { nullptr };

  bool m_makeZClustering {true};  // z_clustering_option
  /* bool do_hit_assoc = true; */
  /* bool do_read_raw = false; */
  /* int m_cluster_version = 4; */
};

#endif

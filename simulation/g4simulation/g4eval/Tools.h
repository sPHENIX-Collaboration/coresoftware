#ifndef CLUSTERMATCHTOOLS__H
#define CLUSTERMATCHTOOLS__H

#include <trackbase/TrkrDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <utility>
#include <array>
#include <set>
/* #include <trackbase_historic/TrackSeed.h> */

class TrkrClusterContainer;
class PHCompositeNode;
class TrkrCluster;
class SvtxTrack;

namespace G4Eval {
  class TrkrClusterComparer {

    // most members are public for easy access after the node has been used
    public:
    TrkrClusterComparer (float _nphi_widths=0.5, float _nz_widths=0.5 );
    int init(PHCompositeNode* topNode, 
        std::string name_truth_clusters="TRKR_TRUTHCLUSTERCONTAINER", 
        std::string name_reco_clusters="TRKR_CLUSTER");

    bool status_good {false};
    TrkrCluster* clus_T {nullptr};
    TrkrCluster* clus_R {nullptr};
    int layer;

    std::pair<bool, float> operator()  // return pair(is_matched,how_good_match)
      (TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R);

    // Members that are set with each set of cluster keys that
    // are passed to it.
    // z and phi locations of phg4 hit (T) and Svtx hit (R)
    float z_T, z_R, phi_T, phi_R; 
    float phisize_R, phisize_T, zsize_R, zsize_T; // sizes are in pixels/layers
    float phi_delta; // abs(phi_T-phi_R)
    float z_delta;   // abs(z_T-z_R)

    bool in_tpc  {false};
    bool in_mvtx {false};
    bool in_intt {false};
    bool in_tpot {false};

    //z pixel sizes. n.b.: there is no z clustering in the INTT
    float m_zstep_tpc  {0.}; // from tpc geometry
    float m_zstep_mvtx {0.};
    // TPOT not implemented yet...

    void set_nz_widths(float   val) { m_nz_widths   = val; };
    void set_nphi_widths(float val) { m_nphi_widths = val; };
    private:
    //phi pixel sizes, got for the geometries from the topNode
    std::array<double, 56> m_phistep {0.}; // the phistep squared
    float m_nphi_widths;
    float m_nz_widths;

    TrkrClusterContainer* m_TruthClusters {nullptr};
    TrkrClusterContainer* m_RecoClusters  {nullptr};
  };
  
  // The following is a struct to iterate over the cluster keys for a given
  // StvxTrack* tracks, starting with the silicone seed and then returning
  // values for the tpc seed. It is used like:
  //
  // for (auto& cluskey : ClusKeyIter(svtx_track)) {
  //    ... // do things with cluster keys
  // }
  struct ClusKeyIter {
    typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
    typedef ClusterKeySet::iterator ClusterKeyIter;

    ClusKeyIter(SvtxTrack* _track);
    // data
    SvtxTrack* track;
    bool in_silicon;
    bool has_tpc;
    bool no_data; // neither a tpc nor a silicon seed
    ClusterKeyIter iter             { };
    ClusterKeyIter iter_end_silicon { };

    ClusKeyIter begin();
    ClusKeyIter end();
    void operator++();
    TrkrDefs::cluskey operator*();
    bool operator!=(const ClusKeyIter& rhs);
  };

  struct HitSetClusKeyVec {
    using Vector = std::vector<std::pair<TrkrDefs::hitsetkey,TrkrDefs::cluskey>>;
    using Iter = Vector::iterator;

    vector keys_svtx;
    vector keys_phg4;

    vector bool matches_svtx;
    vector bool matches_phg4;

    double find_matches(TrkrClusterCommparer&);

    int addSvtxClusters(SvtxTrack*);
    int addPHG4Clusters(TrkrTruthTrack*);

    std::array<unsigned int,4> svtx_cntclus_allMvtxInttTpc();
    std::array<unsigned int,4> phg4_cntclus_allMvtxInttTpc();

    std::array<unsigned int,4> svtx_cntmatchclus_allMvtxInttTpc();
    std::array<unsigned int,4> phg4_cntmatchclus_allMvtxInttTpc();


    // Does the following:
    // for a given SvtxTrack or TrkrTruthTrack,
    // - makes a vector of pairs: {hitsetkey, cluskey}, and sort it
    // - can return pointers to first and last pair with a given hitsetkey
    // - can keep track of which one's are matched and which are not
    HitSetClusKeyVec(SvtxTrack*);
    HitSetClusKeyVec(TrkrTruthTrack*);
    std::vector< std::pair<TrkrDefs::hitsetkey, TrkrDefs::cluskey>> data;
    int cnt(TrkrDefs::hitsetkey); // count how many pairs have 
    vector<TrkrDefs::cluskey

  }

}

#endif

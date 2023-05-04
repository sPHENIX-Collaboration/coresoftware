#ifndef CLUSTERMATCHTOOLS__H
#define CLUSTERMATCHTOOLS__H

#include <trackbase/TrkrDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <utility>
#include <array>
#include <set>
#include <vector>
#include <Eigen/Core>
#include <cfloat>
#include <tuple>
/* #include <trackbase_historic/TrackSeed.h> */

class ActsGeometry;
class EmbRecoMatchContainer;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrTruthTrack;

namespace G4Eval {
  // ClusLoc holds layer, location, phi size and z size
  using ClusLoc = std::tuple<int,Eigen::Vector3d,int,int>;

  // Following function writes msg to the currently active TFile
  // if f_outname is provided, then it will write the message to a new
  // TFiles of that name and close it again.
  void write_StringToTFile(const std::string& msg_name, const std::string& msg);

  std::vector<int> unmatchedSvtxTrkIds(EmbRecoMatchContainer*, SvtxTrackMap*);

  class TrkrClusterComparer {
    // most members are public for easy access after the node has been used
    public:
    TrkrClusterComparer (float _nphi_widths=0.5, float _nz_widths=0.5 );
    int init(PHCompositeNode* topNode, 
        const std::string& name_truth_clusters="TRKR_TRUTHCLUSTERCONTAINER", 
        const std::string&name_reco_clusters="TRKR_CLUSTER");

    TrkrCluster* clus_T { nullptr };
    TrkrCluster* clus_R { nullptr };

    /* std::pair<bool, float> is_match_b */
    std::pair<bool, float> operator() (TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R);

    // Members that are set with each set of cluster keys that
    // are passed to it.
    // z and phi locations of phg4 hit (T) and Svtx hit (R)
    bool is_match  { false };
    int  layer     { INT_MAX };

    float z_T       { FLT_MAX }, z_R       { FLT_MAX };
    float phi_T     { FLT_MAX }, phi_R     { FLT_MAX };
    float phisize_R { FLT_MAX }, phisize_T { FLT_MAX }; // phisize is in nbins * nwidhts
    float zsize_R   { FLT_MAX }, zsize_T   { FLT_MAX }; // zsize   is in nbins * nwdiths
    float phi_delta { FLT_MAX }, z_delta   { FLT_MAX }; // deltas are also in nbins

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

    ClusLoc clusloc_PHG4(std::pair<TrkrDefs::hitsetkey,TrkrDefs::cluskey>);
    ClusLoc clusloc_SVTX(std::pair<TrkrDefs::hitsetkey,TrkrDefs::cluskey>);

    TrkrClusterContainer* m_TruthClusters {nullptr};
    TrkrClusterContainer* m_RecoClusters  {nullptr};
    private:
    //phi pixel sizes, got for the geometries from the topNode
    std::array<double, 56> m_phistep {0.}; // the phistep squared
    float m_nphi_widths;
    float m_nz_widths;

    ActsGeometry*         m_ActsGeometry  {nullptr};

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

  int trklayer_0123(TrkrDefs::hitsetkey); // 0:Mvtx 1:Intt 2:Tpc 3:Tpot 

  class ClusCntr {
    private:
    using Vector = std::vector<std::pair<TrkrDefs::hitsetkey,TrkrDefs::cluskey>>;
    using Iter = Vector::iterator;

    TrkrClusterComparer* comp;
    std::array<int,5> cntclus(Vector& keys);
    std::array<int,5> cnt_matchedclus(Vector& keys, std::vector<bool>& matches);

    public:
    ClusCntr(TrkrClusterComparer*_=nullptr) : comp{_} {};
    TrkrClusterContainer* get_PHG4_clusters();
    TrkrClusterContainer* get_SVTX_clusters();

    Vector svtx_keys {};
    Vector phg4_keys {};

    double match_stat {0};

    void reset();
    std::array<int,3> find_matches();// populated matches_{svtx,phg4}; 
                                     // return's {n-matched, n-phg4, n-svtx}
    std::array<int,3> find_matches(TrkrTruthTrack* g4_track, SvtxTrack* sv_track);

    int phg4_n_matched(); // also same as phg4_cnt_matchedclus()[4]
    int svtx_n_matched(); // should be almost always the same
                     // which is ALMOST guaranteed to be same as svtx_cnt_matchedclus()[4]
    int phg4_nclus() { return (int) phg4_keys.size(); }
    int svtx_nclus() { return (int) svtx_keys.size(); }

    std::vector<bool> svtx_matches;
    std::vector<bool> phg4_matches;

    int addClusKeys(SvtxTrack*); // return number of clusters
    int addClusKeys(TrkrTruthTrack*); // return number of clusters

    std::array<int,5> svtx_cntclus() { return cntclus(svtx_keys); }; // Mvtx Intt Tpc TPOT Sum
    std::array<int,5> phg4_cntclus() { return cntclus(phg4_keys); };

    std::array<int,5> svtx_cnt_matchedclus() {return cnt_matchedclus(svtx_keys, svtx_matches); };
    std::array<int,5> phg4_cnt_matchedclus() {return cnt_matchedclus(phg4_keys, phg4_matches); };

    //I need the cluster widths for diagnostics, too
    std::vector<ClusLoc> phg4_clusloc_all       ();
    std::vector<ClusLoc> phg4_clusloc_unmatched();
    std::vector<ClusLoc> svtx_clusloc_all       ();
    std::vector<ClusLoc> svtx_clusloc_unmatched ();
    std::vector<ClusLoc> clusloc_matched        ();

    void set_comparer(TrkrClusterComparer* _comp) { comp = _comp; };

  };
}

#endif

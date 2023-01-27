#ifndef TRUTHTRKMATCHER__H
#define TRUTHTRKMATCHER__H

#include <fun4all/SubsysReco.h> 
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TTree.h>

#include <array>
#include <iostream>
#include <map>
#include <set>
#include <set>
#include <tuple>
#include <vector>

class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrTruthTrack;
class TrkrTruthTrackContainer;

class TruthRecoTrackMatching : public SubsysReco 
{
  //--------------------------------------------------
  // Standard public interface
  //--------------------------------------------------
  public:
    TruthRecoTrackMatching(                      // Criteria to match a TrkrClusterContainer and track
        /*-----------------------------------------------------------------------------------
         * Input criteria for Truth Track (with nT clusters) to reco track (with nR clusters) :
         *  - nmin_match  : minimum number of clusters in Truth and Reco track that must 
         *                  match for tracks to match
         *  - cutoff_dphi : maximum distance out in |phi_truth-phi_reco| for a matched track
         *  - same_dphi   : within this distance, all tracks must be checked
         *  - cutoff_deta :   like cutoff_dphi but for eta
         *  - same_deta   :   like same_dphi but for eta
         *  - cluster_nzwidths   : z-distance allowed for the truth center to be from 
         *                         reco-cluster center:
         *                         |z_true-z_reco| must be <= dz_reco * cluster_nzwidths
         *  - cluster_nphiwidths : like cluster_nzwidths for phi
         *--------------------------------------------------------*/
          const unsigned short _nmin_match = 4    
        , const float  _nmin_ratio         = 0.   
        , const double _cutoff_dphi        = 0.3  
        , const double _same_dphi          = 0.05 
        , const double _cutoff_deta        = 0.3  
        , const double _same_deta          = 0.05 
        , const double _cluster_nzwidths   = 0.5
        , const double _cluster_nphiwidths = 0.5
        , const unsigned short    _max_nreco_per_truth  = 1
        , const unsigned short    _max_ntruth_per_reco  = 1
    );  // for some output kinematis

    ~TruthRecoTrackMatching() override = default; 
    int Init(PHCompositeNode          *) override  { return 0; };
    int InitRun(PHCompositeNode       *) override; //`
    int process_event(PHCompositeNode *) override; //`
    int End(PHCompositeNode           *) override;

    int createNodes(PHCompositeNode* topNode);

    void set_cluster_nphiwidths       (float val) { m_cluster_nphiwidths = val; };
    void set_cluster_nzwidths         (float val) { m_cluster_nzwidths   = val; };
    void set_cutoff_deta              (float val) { m_cutoff_deta        = val; };
    void set_cutoff_dphi              (float val) { m_cutoff_dphi        = val; };
    void set_nmin_truth_cluster_ratio (float val) { m_nmincluster_ratio  = val; };
    void set_smallsearch_deta         (float val) { m_same_deta          = val; };
    void set_smallsearch_dphi         (float val) { m_same_dphi          = val; };

    void set_max_nreco_per_truth (unsigned short val) { m_max_nreco_per_truth = val; };
    void set_max_ntruth_per_reco (unsigned short val) { m_max_ntruth_per_reco = val; };

  private:
  //--------------------------------------------------
  // Internal functions
  //--------------------------------------------------

    //--------------------------------------------------
    // Constant parameters for track matching
    //--------------------------------------------------
    unsigned short  m_nmincluster_match; // minimum of matched clustered to keep a truth to emb match
    float  m_nmincluster_ratio; // minimum ratio of truth clustered that must be matched in reconstructed track

    double m_cutoff_dphi; // how far in |phi_truth-phi_reco| to match
    double m_same_dphi;   // |phi_truth-phi_reco| to auto-evaluate (is < _m_cutoff_dphi)
    double m_cutoff_deta; //  how far in |eta_truth-eta_reco| to match
    double m_same_deta;   // |eta_truth-eta_reco| to auto-evaluate (is < m_cutoff_deta)

    double m_cluster_nzwidths;   // cutoff in *getPhiSize() in cluster for |cluster_phi_truth-cluster_phi_reco| to match
    double m_cluster_nphiwidths; // same for eta

    unsigned short m_max_nreco_per_truth;
    unsigned short m_max_ntruth_per_reco;


    std::array<double, 55> m_phistep {0.}; // the phistep squared
    double m_zstep {0.};

    std::map<unsigned short, unsigned short>  m_nmatched_index_true {};

    std::map<unsigned short, unsigned short>* m_nmatched_id_reco {nullptr};
    std::map<unsigned short, unsigned short>* m_nmatched_id_true {nullptr};


    //--------------------------------------------------
    // Data from input nodes
    //--------------------------------------------------
    PHG4TruthInfoContainer       *m_PHG4TruthInfoContainer       {nullptr}; // Get the truth track ids
    SvtxTrackMap                 *m_SvtxTrackMap                 {nullptr};
    TrkrClusterContainer         *m_TruthClusterContainer        {nullptr};
    TrkrClusterContainer         *m_RecoClusterContainer         {nullptr};
    TrkrTruthTrackContainer      *m_TrkrTruthTrackContainer      {nullptr};
    PHG4TpcCylinderGeomContainer *m_PHG4TpcCylinderGeomContainer {nullptr};

    ActsGeometry* m_ActsGeometry {nullptr};

    // Output data node:
    EmbRecoMatchContainer   *m_EmbRecoMatchContainer   {nullptr};

    //--------------------------------------------------
    //    RECO data for a "table" of reconstructed tracks
    //--------------------------------------------------
    using RECOentry = std::tuple<float, float, float, unsigned short>;
    static constexpr int RECOphi = 0;
    static constexpr int RECOeta = 1;
    static constexpr int RECOpt  = 2;
    static constexpr int RECOid  = 3;
  public:
    using RECOvec   = std::vector<RECOentry>;
    using RECOiter = RECOvec::iterator;
    using RECO_pair_iter = std::pair<RECOiter, RECOiter>;
  private:

    struct CompRECOtoPhi {
      bool operator() (const RECOentry& lhs, const double&    rhs) { return std::get<RECOphi>(lhs) < rhs; }
      bool operator() (const double&    lhs, const RECOentry& rhs) { return lhs < std::get<RECOphi>(rhs); }
      bool operator() (const RECOentry& lhs, const RECOentry& rhs) { return std::get<RECOphi>(lhs) < std::get<RECOphi>(rhs); }
    };
    struct CompRECOtoEta {
      bool operator() (const RECOentry& lhs, const double&    rhs) { return std::get<RECOeta>(lhs) < rhs; }
      bool operator() (const double&    lhs, const RECOentry& rhs) { return lhs < std::get<RECOeta>(rhs); }
      bool operator() (const RECOentry& lhs, const RECOentry& rhs) { return std::get<RECOeta>(lhs) < std::get<RECOeta>(rhs); }
    };
    struct CompRECOtoPt  {
      bool operator() (const RECOentry& lhs, const double&    rhs) { return std::get<RECOpt>(lhs) < rhs; }
      bool operator() (const double&    lhs, const RECOentry& rhs) { return lhs < std::get<RECOpt>(rhs); }
      bool operator() (const RECOentry& lhs, const RECOentry& rhs) { return std::get<RECOpt>(lhs) < std::get<RECOpt>(rhs); }
    };

    // sorting tuple is by: phi, eta, pT,  index, is_matched
    //--------------------------------------------------
    //    PossibleMatches (just array<unsinged short, 5>)
    //--------------------------------------------------
  public:
    using PossibleMatch = std::array<unsigned short, 5>;
  private:
    static constexpr int PM_nmatch  = 0;
    static constexpr int PM_ntrue  = 1;
    static constexpr int PM_nreco  = 2;
    static constexpr int PM_idtrue = 3;
    static constexpr int PM_idreco = 4;
    struct SortPossibleMatch  {
      // Sort by most matched clusters first, then smallest number of truth clusters, then smallest number of reco clusters
      bool operator() (const PossibleMatch& lhs, const PossibleMatch& rhs) { 
        if (lhs[PM_nmatch] != rhs[PM_nmatch]) return lhs[PM_nmatch] > rhs[PM_nmatch];
        if (lhs[PM_ntrue ] != rhs[PM_ntrue ]) return lhs[PM_ntrue ] < rhs[PM_ntrue ];
        if (lhs[PM_nreco ] != rhs[PM_nreco ]) return lhs[PM_nreco ] < rhs[PM_nreco ];
        return false;
      }
    };

    //--------------------------------------------------
    //    PossibleMatches (just array<unsinged short, 5>)
    //--------------------------------------------------
  public:
    //--------------------------------------------------
    // non-node member data
    //--------------------------------------------------
    RECOvec recoData {}; // "sv" = Sort Vector

    //--------------------------------------------------
    // Member functions
    //--------------------------------------------------
   
    //    Main functions
    // -------------------------------------------------------------------
    std::pair<std::vector<unsigned short>, std::vector<unsigned short>> 
      find_box_matches(float truth_phi, float truth_eta, float truth_pt); // will populate to truth_to_reco_map and 
    void match_tracks_in_box( std::vector<std::pair<unsigned short,unsigned short>>& indices );  // pairs of {id_true, id_reco}
    
    //    Helper functions
    // -------------------------------------------------------------------
    float delta_outer_pt (float) const;
    float delta_inner_pt (float) const;
    float abs_dphi (float phi0, float phi1);
    float sigma_CompMatchClusters (PossibleMatch&);
    bool  skip_match (PossibleMatch& match);
    bool  at_nmax_index_true (unsigned short) ; // test if the truth track already has maximum matches matches
    bool  at_nmax_id_reco    (unsigned short)    ;    // "           reco                                          "

    std::pair<bool, float> compare_cluster_pair(TrkrDefs::cluskey key_T,
        TrkrDefs::cluskey key_R, TrkrDefs::hitsetkey key, bool
        calc_sigma=false);
};



#endif

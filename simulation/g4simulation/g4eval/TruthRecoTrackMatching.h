#ifndef TRUTHTRKMATCHER__H
#define TRUTHTRKMATCHER__H

#include <trackbase/TrkrDefs.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <tuple>
#include <vector>
#include <map>
#include <array>
#include <set>
#include <iostream>
#include <TTree.h>
#include <TFile.h>

#include <trackbase/ActsGeometry.h>
/* #include <Acts/Definitions/Algebra.hpp> */
/* #include <Eigen/Core> */
/* #include <Eigen/Dense> */


class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxTrack;
class TrackSeedContainer;
class TrkrClusterContainer;
class TrkrCluster;
class TrkrTruthTrackContainer;
class TrkrTruthTrack;

class PHG4TpcCylinderGeom;
class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;

class TruthRecoTrackMatching : public SubsysReco 
{
  //--------------------------------------------------
  // Standard public interface
  //--------------------------------------------------
  public:
    TruthRecoTrackMatching(                      // Criteria to match a TrkrClusterContainer and track
        const unsigned short _nmin_match = 4,    // Min. matched clusters to match matchs
        const float  _nmin_ratio         = 0.,   // "                  " ratio to clusters in truth to match tracks
        const double _cutoff_dphi        = 0.3,  // |phi_true-phi_reco| bounds to match tracks (bigger)
        const double _same_dphi          = 0.05, // |phi_true-phi_reco| bounds *must* check track
        const double _cutoff_deta        = 0.3,  //  same for eta as for phi
        const double _same_deta          = 0.05, //  '                     '
        const double _cluster_nzwidths   = 1.,   //  to match clusters: max |z_true-z_reco| in ratio to sqrt(dz_t^2+dz_m^2)
        const double _cluster_nphiwidths = 1.,
        const std::string _m_filename = "" );  // same for phi

    ~TruthRecoTrackMatching() override = default; 
    int Init(PHCompositeNode          *) override { return 0; };
    int InitRun(PHCompositeNode       *) override; //`
    int process_event(PHCompositeNode *) override; //`
    int End(PHCompositeNode           *) override;

  private:
  //--------------------------------------------------
  // Internal functions
  //--------------------------------------------------

    //--------------------------------------------------
    // Constant parameters for track matching
    //--------------------------------------------------
    const unsigned short  m_nmincluster_match; // minimum of matched clustered to keep a truth to emb match
    const float  m_nmincluster_ratio; // minimum ratio of truth clustered that must be matched in reconstructed track

    const double m_cutoff_dphi; // how far in |phi_truth-phi_reco| to match
    const double m_same_dphi;   // |phi_truth-phi_reco| to auto-evaluate (is < _m_cutoff_dphi)
    const double m_cutoff_deta; //  how far in |eta_truth-eta_reco| to match
    const double m_same_deta;   // |eta_truth-eta_reco| to auto-evaluate (is < m_cutoff_deta)

    const double m_cluster_nzwidths_2;   // cutoff in *getPhiSize() in cluster for |cluster_phi_truth-cluster_phi_reco| to match
    const double m_cluster_nphiwidths_2; // same for eta

    std::array<double, 55> m_phistep_2; // the phistep squared
    double m_zstep_2;

    //--------------------------------------------------
    // Data from input nodes
    //--------------------------------------------------
    PHG4TruthInfoContainer       *m_PHG4TruthInfoContainer  {nullptr}; // Get the truth track ids
    SvtxTrackMap                 *m_SvtxTrackMap            {nullptr};
    TrkrClusterContainer         *m_TruthClusterContainer   {nullptr};
    TrkrClusterContainer         *m_RecoClusterContainer    {nullptr};
    TrkrTruthTrackContainer      *m_TrkrTruthTrackContainer {nullptr};
    PHG4TpcCylinderGeomContainer *m_PHG4TpcCylinderGeomContainer  {nullptr};

    ActsGeometry* m_ActsGeometry {nullptr};

    // Output data node:
    EmbRecoMatchContainer   *m_EmbRecoMatchContainer   {nullptr};

    // Output TTree for the sake of checking results statistics:
    TFile * m_fileout;
    TTree * m_treeout;
    bool    m_maketree { true };

    // branches for the output tree
    float b_reco_phi,  b_reco_eta,  b_reco_pt,
          b_truth_phi, b_truth_eta, b_truth_pt;

    std::vector<int>   b_reco_cl_layer, b_truth_cl_layer;

    std::vector<float> b_reco_cl_phi,  b_reco_cl_phiwidth,  b_reco_cl_z,  b_reco_cl_zwidth,
                       b_truth_cl_phi, b_truth_cl_phiwidth, b_truth_cl_z, b_truth_cl_zwidth;

    //----------------------------------------------------------
    // Local data types for manipulation and passing
    //   * provide indices to arrays and tuples in constexp ints
    //----------------------------------------------------------

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

    // sorting tuple is by: is_matched, phi, eta, pT,  index, is_matched
    // c++ will be default sort these by index in order

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
    using CompCluster = std::array<float,6>;
  private:
    static constexpr int COMP_dphi = 0;
    static constexpr int COMP_phisize_true = 1; // note this is converted with step size already
    static constexpr int COMP_phisize_reco = 2; // converted with stepsize
    static constexpr int COMP_dz = 3;
    static constexpr int COMP_zsize_true = 4; // converted with stepsize
    static constexpr int COMP_zsize_reco = 5; // converted with stepsize
  
    //--------------------------------------------------
    // non-node member data
    //--------------------------------------------------
    RECOvec recoData {}; // "sv" = Sort Vector

    //--------------------------------------------------
    // Member functions
    //--------------------------------------------------
   
    //    Main functions
    // -------------------------------------------------------------------
    int createNodes(PHCompositeNode* topNode);
    std::pair<std::vector<unsigned short>, std::vector<unsigned short>> 
      find_box_matches(float truth_phi, float truth_eta, float truth_pt); // will populate to truth_to_reco_map and 
    unsigned int match_tracks_in_box( std::vector<std::pair<unsigned short,unsigned short>>& indices,  // pairs of {id_true, id_reco}
      std::set<int>& true_matched,
      std::set<int>& reco_matched,
      bool update_matched_sets);
    
    //    Helper functions
    // -------------------------------------------------------------------
    float delta_outer_pt(float) const;
    float delta_inner_pt(float) const;
    float abs_dphi (float phi0, float phi1);

    //       Functions for gettings sets of TrkrClusters associated with tracks
    // ------------------------------------------------------------------------
    std::array<TrkrDefs::cluskey,55> truekey_arr55   (int id_true);              // get cluster keys for truth track
    std::array<std::vector<TrkrDefs::cluskey>,55> recokey_arr55vec(int id_reco); // get cluster keys for reco track

    PossibleMatch make_PossibleMatch(
        unsigned short index_reco, 
        unsigned short index_truth, 
        std::array<TrkrDefs::cluskey,55>& keys_truth);

    /* CompCluster                         make_CompCluster(TrkrDefs::cluskey truthkey, TrkrDefs::cluskey recokey); */
    float sigma_CompCluster(PossibleMatch&);
    bool skip_match(PossibleMatch& match, std::set<int>& true_matched, std::set<int>& reco_matched); 
    // skip match is true or reco track already matched
    std::pair<bool, double> compare_clusters(TrkrDefs::cluskey truthkey, TrkrDefs::cluskey recokey, bool print=false);
    /* bool is_match(const CompCluster&); */
    bool match_pass_cuts(PossibleMatch& match);
  
    // locally - maybe?
    /* TrackSeedContainer *m_tpcTrackSeedContainer   {nullptr}; // Get the seeds from the tracks to get the clusters */
};



#endif

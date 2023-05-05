#ifndef FILLTRUTHRECOMATCHTREE_H
#define FILLTRUTHRECOMATCHTREE_H

/**
 * @file trackbase/TrkrMatchDefs.h
 * @author D. Stewart
 * @date February 2023
 * @brief Write a TFile with a TTree with the matched tracks and (optionall) clusters
 *   This module takes the matching done in TruthRecoTrackMatching in the 
 *   EmbRecoMatchContainer objects, and writes out the tracks kinematics to TTree
 *   in a TFile.
 *    The following data is always written out:
 *      For each track type:
 *          G4M : PHG4 matched
 *          G4U : PHG4 unmatched
 *          SvM : Svtx (reco) matched
 *          SvU : SVtx (reco) unmatched
 *      The following is stored:
 *          trackid, nclus, nclus{mvtx,intt,tpc}, pt, phi, eta
 *      For matched tracks (G4M and SvM):
 *           nclus_matchrat, nclus{mvtx,intt,tpc}_matchrat
 *      For matched tracks, the numbers of matched clusters (these are shared by G4M and SvM):
 *           nclusM, nclusM_{mvtx,intt,tpc}
 *      
 *    If the option is passed to save clusteres (default is on), then the cluster
 *    locations are saved too, as unmatched clusters (in all types of tracks),
 *    and the mutually matched clusters (shared by the G4M and SvM tracks):
 *      {G4U,G4M,SvU,SvM}_clusU_{i0,i1,x,y,z,r}
 *    For matche clusters, these are simply:
 *      clusM_{i0,i1,x,y,z,r,layer}
 *    The vectors of x,y,z,r are stoped in each branch with each tracks' data sequentially
 *    following the last. The branch i0 and i1 index where the one starts and the other stops.
 *    for each track.
 *
 *    Data for matched tracks positionally align with each other (i.e. 1st entry in each
 *    correlate, then the 2nd, etc...)
 *
 *    A track is only an unmatched track is it has no matches (the options in 
 *    TruthRecoTrackMatching can allow multiply matches for the same track...)
 *
 *    options:
 *      save cluster locations = true
 *      save un-matched Svtx tracks = false
 */

#include "g4evaltools.h"

#include <fun4all/SubsysReco.h> 
#include <vector>
#include <string>

class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4ParticleSvtxMap;
class SvtxPHG4ParticleMap;
class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeedContainer;
class TrkrCluster;
/* class TrkrClusterContainer; */
class TrkrTruthTrack;
class TrkrTruthTrackContainer;
class TTree;
class TH2D;

class FillTruthRecoMatchTree : public SubsysReco
{
 public:
  FillTruthRecoMatchTree(
        bool  _fill_clusters      = true
      , bool  _fill_svtxnomatch   = false
      , float _cluster_nzwidths   = 0.5
      , float _cluster_nphiwidths = 0.5
      , const std::string& tfile_name="trackclusmatch.root"
  );

  virtual ~FillTruthRecoMatchTree();

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override;
  int End(PHCompositeNode *topNode) override;

  void clear_clusvecs(std::string tag="");

  void print_mvtx_diagnostics();

 private:

   int createNodes(PHCompositeNode *topNode);

   G4Eval::TrkrClusterComparer m_cluster_comp { 1., 1.};
   G4Eval::ClusCntr            m_cluscntr;

   // contianer used to fill the other track matches
   EmbRecoMatchContainer        *m_EmbRecoMatchContainer   {nullptr}; 
   PHG4TruthInfoContainer       *m_PHG4TruthInfoContainer       {nullptr};
   SvtxTrackMap                 *m_SvtxTrackMap                 {nullptr};
   /* TrkrClusterContainer         *m_TruthClusterContainer        {nullptr}; */
   /* TrkrClusterContainer         *m_RecoClusterContainer         {nullptr}; */
   TrkrTruthTrackContainer      *m_TrkrTruthTrackContainer      {nullptr};
   /* PHG4TpcCylinderGeomContainer *m_PHG4TpcCylinderGeomContainer {nullptr}; */


   TTree* m_ttree;
   bool   m_fill_clusters;
   bool   m_fill_SvU; // unmatched Svtx tracks
   std::string m_outfile_name;

  // Tree Branch members:
    int   nevent       {-1};
    int   nphg4        {0};
    int   nsvtx        {0};
    int   ntrackmatches {0};
    int   nphg4_part   {0};
    float centrality   {0.};

    // Tracks and clustes
    //
    // lables:
    //  tracks:
    //    g4 : phg4track matched
    //    sv : svtx_track matched
    //    gU : phg4track not-matched
    //    sU : svtx_track not-matched
    //  clusters:
    //    M : matched
    //    U : unmatched
    TH2D* h2_G4_nPixelsPhi;
    TH2D* h2_G4_nPixelsZ;
    TH2D* h2_Sv_nPixelsPhi;
    TH2D* h2_Sv_nPixelsZ;
    
    // Track tree
    int  b_trackid;
    bool b_is_g4track;
    bool b_is_Svtrack;
    bool b_is_matched;

    float b_trkpt;
    float b_trkphi;
    float b_trketa;


    int   b_nclus         {};
    int   b_nclustpc      {};
    int   b_nclusmvtx     {};
    int   b_nclusintt     {};

    float b_matchrat      {};
    float b_matchrat_intt {};
    float b_matchrat_mvtx {};
    float b_matchrat_tpc  {};

    std::vector<bool>  b_clusmatch  {};
    std::vector<float> b_clus_x     {};
    std::vector<float> b_clus_y     {};
    std::vector<float> b_clus_z     {};
    std::vector<float> b_clus_r     {};
    std::vector<int>   b_clus_layer {};
    std::vector<int>   b_clus_nphibins   {};
    std::vector<int>   b_clus_ntbins     {};

    /* std::vector<int>   b_G4M_trackid            {}; // g4-track-matched */
    /* std::vector<int>   b_G4M_nclus              {}; */
    /* std::vector<int>   b_G4M_nclusmvtx          {}; */
    /* std::vector<int>   b_G4M_nclusintt          {}; */
    /* std::vector<int>   b_G4M_nclustpc           {}; */
    /* std::vector<float> b_G4M_nclus_matchrat     {}; */
    /* std::vector<float> b_G4M_nclusmvtx_matchrat {}; */
    /* std::vector<float> b_G4M_nclusintt_matchrat {}; */
    /* std::vector<float> b_G4M_nclustpc_matchrat  {}; */
    /* std::vector<float> b_G4M_pt                 {}; */
    /* std::vector<float> b_G4M_phi                {}; */
    /* std::vector<float> b_G4M_eta                {}; */
    /* std::vector<int>   b_SvM_trackid            {}; // Svtx-track-matched */
    /* std::vector<int>   b_SvM_nclus              {}; */
    /* std::vector<int>   b_SvM_nclusmvtx          {}; */
    /* std::vector<int>   b_SvM_nclusintt          {}; */
    /* std::vector<int>   b_SvM_nclustpc           {}; */
    /* std::vector<float> b_SvM_nclus_matchrat     {}; */
    /* std::vector<float> b_SvM_nclusmvtx_matchrat {}; */
    /* std::vector<float> b_SvM_nclusintt_matchrat {}; */
    /* std::vector<float> b_SvM_nclustpc_matchrat  {}; */
    /* std::vector<float> b_SvM_pt                 {}; */
    /* std::vector<float> b_SvM_phi                {}; */
    /* std::vector<float> b_SvM_eta                {}; */
    /* std::vector<int>   b_clusM_i0               {}; // if storing clusters -- matched clusters */
    /* std::vector<int>   b_clusM_i1               {}; */
    /* std::vector<float> b_clusM_layer            {}; */
    /* std::vector<float> b_clusM_x                {}; */
    /* std::vector<float> b_clusM_y                {}; */
    /* std::vector<float> b_clusM_z                {}; */
    /* std::vector<float> b_clusM_r                {}; */
    /* std::vector<int>   b_G4M_clusU_i0           {}; // matched phg4 unmatched clusters */
    /* std::vector<int>   b_G4M_clusU_i1           {}; */
    /* std::vector<float> b_G4M_clusU_layer        {}; */
    /* std::vector<float> b_G4M_clusU_x            {}; */
    /* std::vector<float> b_G4M_clusU_y            {}; */
    /* std::vector<float> b_G4M_clusU_z            {}; */
    /* std::vector<float> b_G4M_clusU_r            {}; */
    /* std::vector<int>   b_SvM_clusU_i0           {}; //  matched phg4 unmatched clusters */
    /* std::vector<int>   b_SvM_clusU_i1           {}; */
    /* std::vector<float> b_SvM_clusU_layer        {}; */
    /* std::vector<float> b_SvM_clusU_x            {}; */
    /* std::vector<float> b_SvM_clusU_y            {}; */
    /* std::vector<float> b_SvM_clusU_z            {}; */
    /* std::vector<float> b_SvM_clusU_r            {}; */
    /* std::vector<int>   b_G4U_trackid            {}; // unmatched tracks */
    /* std::vector<int>   b_G4U_nclus              {}; */
    /* std::vector<int>   b_G4U_nclusmvtx          {}; */
    /* std::vector<int>   b_G4U_nclusintt          {}; */
    /* std::vector<int>   b_G4U_nclustpc           {}; */
    /* std::vector<float> b_G4U_pt                 {}; */
    /* std::vector<float> b_G4U_phi                {}; */
    /* std::vector<float> b_G4U_eta                {}; */
    /* std::vector<int>   b_SvU_trackid            {}; // Svtx-track-matched */
    /* std::vector<int>   b_SvU_nclus              {}; */
    /* std::vector<int>   b_SvU_nclusmvtx          {}; */
    /* std::vector<int>   b_SvU_nclusintt          {}; */
    /* std::vector<int>   b_SvU_nclustpc           {}; */
    /* std::vector<float> b_SvU_pt                 {}; */
    /* std::vector<float> b_SvU_phi                {}; */
    /* std::vector<float> b_SvU_eta                {}; */
    /* std::vector<int>   b_G4U_clusU_i0           {}; // unmatched phg4 unmatched clusters */
    /* std::vector<int>   b_G4U_clusU_i1           {}; */
    /* std::vector<float> b_G4U_clusU_layer        {}; */
    /* std::vector<float> b_G4U_clusU_x            {}; */
    /* std::vector<float> b_G4U_clusU_y            {}; */
    /* std::vector<float> b_G4U_clusU_z            {}; */
    /* std::vector<float> b_G4U_clusU_r            {}; */
    /* std::vector<int>   b_SvU_clusU_i0           {}; // unmatched phg4 unmatched clusters */
    /* std::vector<int>   b_SvU_clusU_i1           {}; */
    /* std::vector<float> b_SvU_clusU_layer        {}; */
    /* std::vector<float> b_SvU_clusU_x            {}; */
    /* std::vector<float> b_SvU_clusU_y            {}; */
    /* std::vector<float> b_SvU_clusU_z            {}; */
    /* std::vector<float> b_SvU_clusU_r            {}; */



};

#endif  // FILLTRUTHRECOMATCHTREE_H

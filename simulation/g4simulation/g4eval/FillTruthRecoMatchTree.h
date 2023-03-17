#ifndef FILLTRUTHRECOMATCHTREE_H
#define FILLTRUTHRECOMATCHTREE_H

/**
 * @file trackbase/TrkrMatchDefs.h
 * @author D. Stewart
 * @date February 2023
 * @brief Write a TFile with a TTree with the matched tracks and (optionall) clusters
 *
 * read out the matched tracks from the node tree,
 * look at the matched (and un-matched) tracks and plot out their
 * kinematics and the locations (potentially) of the matched clusters.
 *
 *    An entry in SvtxPHG4ParticleMap could be:
 *      12 -> 41.46 -> { 2, 4 }
 *         -> 18.46 -> { 7 }
 *    which is to say, reco track id 12 has matches to truth tracks 2, 4, and
 *    7. Matches 12->2 and 12->4 weighting key (41.46) indicate that there were
 *    41 matched clusters, that the reco track had 46 clusters. Match 12->7 key
 *    (18.46) indicates that there were 18 matched clusters (and, again, the reco track had 46 clusters)
 *
 *    Assuming that truth tracks 2, 4, and 7, were matched only to reco track 12, and 
 *    each had 45, 44, and 47 clusters, respectively, then the corresponding entries
 *    in PHG4ParticleSvtxMap would be:
 *      2 -> 41.45 { 12 } 
 *      4 -> 41.44 { 12 } 
 *      7 -> 18.47 { 12 }
 *
 */

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
class TrkrClusterContainer;
class TrkrTruthTrack;
class TrkrTruthTrackContainer;
class TTree;
class TFile;

class FillTruthRecoMatchTree : public SubsysReco
{
 public:
  FillTruthRecoMatchTree(const std::string &name = "FillTruthRecoMatchTree",
      bool _fill_clusters = false, const std::string tfile_name="trackclusmatch.root");

  virtual ~FillTruthRecoMatchTree();

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override;
  int End(PHCompositeNode *topNode) override;

  void clear_vectors();

 private:

   int createNodes(PHCompositeNode *topNode);
    // if looking for matched and un-matched clusters, need
    // the matching criteria.
    // This should match that in TruthRecoTrackMatching
   double _cluster_nzwidths   { 0.5 };
   double _cluster_nphiwidths { 0.5 };

   
   // contianer used to fill the other track matches
   EmbRecoMatchContainer   *m_EmbRecoMatchContainer   {nullptr}; 
   PHG4TruthInfoContainer       *m_PHG4TruthInfoContainer       {nullptr};
   SvtxTrackMap                 *m_SvtxTrackMap                 {nullptr};
   TrkrClusterContainer         *m_TruthClusterContainer        {nullptr};
   TrkrClusterContainer         *m_RecoClusterContainer         {nullptr};
   TrkrTruthTrackContainer      *m_TrkrTruthTrackContainer      {nullptr};
   PHG4TpcCylinderGeomContainer *m_PHG4TpcCylinderGeomContainer {nullptr};

   TFile* m_tfile;
   TTree* m_ttree;
   bool   m_fill_clusters { false };

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
    
    // TRACKS WHICH ARE MATCHED
    std::vector<int>   b_g4M_trackid            {}; // g4-track-matched
    std::vector<int>   b_g4M_nclus              {};
    std::vector<int>   b_g4M_nclusmvtx          {};
    std::vector<int>   b_g4M_nclusintt          {};
    std::vector<int>   b_g4M_nclustpc           {};
    std::vector<float> b_g4M_nclusmvtx_matchrat {};
    std::vector<float> b_g4M_nclusintt_matchrat {};
    std::vector<float> b_g4M_nclustpc_matchrat  {};
    std::vector<float> b_g4M_pt                 {};
    std::vector<float> b_g4M_px                 {};
    std::vector<float> b_g4M_py                 {};
    std::vector<float> b_g4M_pz                 {};
    std::vector<int>   b_g4M_clusM_i0           {}; // g4-track-cluster-matched
    std::vector<int>   b_g4M_clusM_i1           {};
    std::vector<float> b_g4M_clusM_layer        {};
    std::vector<float> b_g4M_clusM_x            {};
    std::vector<float> b_g4M_clusM_y            {};
    std::vector<float> b_g4M_clusM_z            {};
    std::vector<float> b_g4M_clusM_r            {};
    std::vector<int>   b_g4M_clusU_i0           {}; // g4-track-cluster-unmatched
    std::vector<int>   b_g4M_clusU_i1           {};
    std::vector<float> b_g4M_clusU_layer        {};
    std::vector<float> b_g4M_clusU_x            {};
    std::vector<float> b_g4M_clusU_y            {};
    std::vector<float> b_g4M_clusU_z            {};
    std::vector<float> b_g4M_clusU_r            {};
    std::vector<int>   b_svM_trackid            {}; //svtx-track-matched
    std::vector<int>   b_svM_nclus              {};
    std::vector<int>   b_svM_nclusmvtx          {};
    std::vector<int>   b_svM_nclusintt          {};
    std::vector<int>   b_svM_nclustpc           {};
    std::vector<float> b_svM_nclusmvtx_matchrat {};
    std::vector<float> b_svM_nclusintt_matchrat {};
    std::vector<float> b_svM_nclustpc_matchrat  {};
    std::vector<float> b_svM_pt                 {};
    std::vector<float> b_svM_px                 {};
    std::vector<float> b_svM_py                 {};
    std::vector<float> b_svM_pz                 {};
    std::vector<int>   b_svM_clusM_i0           {}; // svtx-track-matched-cluster-matched
    std::vector<int>   b_svM_clusM_i1           {};
    std::vector<float> b_svM_clusM_layer        {};
    std::vector<float> b_svM_clusM_x            {};
    std::vector<float> b_svM_clusM_y            {};
    std::vector<float> b_svM_clusM_z            {};
    std::vector<float> b_svM_clusM_r            {};
    std::vector<int>   b_svM_clusU_i0           {};
    std::vector<int>   b_svM_clusU_i1           {};
    std::vector<float> b_svM_clusU_layer        {};
    std::vector<float> b_svM_clusU_x            {};
    std::vector<float> b_svM_clusU_y            {};
    std::vector<float> b_svM_clusU_z            {};
    std::vector<float> b_svM_clusU_r            {};
    std::vector<int>   b_nclus_match            {}; // Ratios for clusters for matched tracks
    std::vector<int>   b_g4svmatch_nclusmvtx    {};
    std::vector<int>   b_g4svmatch_nclusintt    {};
    std::vector<int>   b_g4svmatch_nclustpc     {};
    std::vector<int>   b_g4U_trackid            {}; // g4-track-unmatched-
    std::vector<int>   b_g4U_nclus              {};
    std::vector<int>   b_g4U_nclusmvtx          {};
    std::vector<int>   b_g4U_nclusintt          {};
    std::vector<int>   b_g4U_nclustpc           {};
    std::vector<float> b_g4U_pt                 {};
    std::vector<float> b_g4U_px                 {};
    std::vector<float> b_g4U_py                 {};
    std::vector<float> b_g4U_pz                 {};
    std::vector<int>   b_g4U_clusU_i0           {}; // g4-track-unmatched-clust-unmatched
    std::vector<int>   b_g4U_clusU_i1           {};
    std::vector<float> b_g4U_clusU_layer        {};
    std::vector<float> b_g4U_clusU_x            {};
    std::vector<float> b_g4U_clusU_y            {};
    std::vector<float> b_g4U_clusU_z            {};
    std::vector<float> b_g4U_clusU_r            {};
    std::vector<int>   b_svU_trackid            {}; // svtx-track-unmatched
    std::vector<int>   b_svU_nclus              {};
    std::vector<int>   b_svU_nclusmvtx          {};
    std::vector<int>   b_svU_nclusintt          {};
    std::vector<int>   b_svU_nclustpc           {};
    std::vector<float> b_svU_pt                 {};
    std::vector<float> b_svU_px                 {};
    std::vector<float> b_svU_py                 {};
    std::vector<float> b_svU_pz                 {};
    std::vector<int>   b_svU_clusU_i0           {}; // svtx-track-unmatched-cluster-unmatched
    std::vector<int>   b_svU_clusU_i1           {};
    std::vector<float> b_svU_clusU_layer        {};
    std::vector<float> b_svU_clusU_x            {};
    std::vector<float> b_svU_clusU_y            {};
    std::vector<float> b_svU_clusU_z            {};
    std::vector<float> b_svU_clusU_r            {};
};

#endif  // FILLTRUTHRECOMATCHTREE_H

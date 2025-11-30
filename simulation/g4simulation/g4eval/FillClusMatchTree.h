#ifndef G4EVAL_FILLCLUSMATCHTREE_H
#define G4EVAL_FILLCLUSMATCHTREE_H

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

#include "TrackClusEvaluator.h"
#include "TrkrClusLoc.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4ParticleSvtxMap;
class SvtxPHG4ParticleMap;
class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4TpcGeom;
class PHG4TpcGeomContainer;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeedContainer;
class TrkrCluster;
class TrkrTruthTrack;
class TrkrTruthTrackContainer;
class TTree;
class TH2;
class ClusHitsVerbose;
class TrkrClusterIsMatcher;

class FillClusMatchTree : public SubsysReco
{
 public:
  FillClusMatchTree(
      TrkrClusterIsMatcher *_ismatcher, const std::string &tfile_name = "trackclusmatch.root", bool _fill_clusters = true, bool _fill_clusverbose = true, bool _fill_SvUnMatched = false);

  virtual ~FillClusMatchTree() = default;

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override;
  int End(PHCompositeNode *topNode) override;

  void clear_clusvecs(const std::string &tag = "");

  void print_mvtx_diagnostics();

  TrkrClusterIsMatcher *m_ismatcher;

 private:
  int createNodes(PHCompositeNode *topNode);
  TrackClusEvaluator m_TCEval;

  // contianer used to fill the other track matches
  EmbRecoMatchContainer *m_EmbRecoMatchContainer{nullptr};
  PHG4TruthInfoContainer *m_PHG4TruthInfoContainer{nullptr};
  SvtxTrackMap *m_SvtxTrackMap{nullptr};
  TrkrTruthTrackContainer *m_TrkrTruthTrackContainer{nullptr};
  ClusHitsVerbose *m_PHG4ClusHitVerb{nullptr};
  ClusHitsVerbose *m_SvtxClusHitVerb{nullptr};
  /* TrkrClusterContainer         *m_TruthClusterContainer        {nullptr}; */
  /* TrkrClusterContainer         *m_RecoClusterContainer         {nullptr}; */
  /* PHG4TpcGeomContainer *m_PHG4TpcGeomContainer {nullptr}; */

  TTree *m_ttree;
  std::string m_outfile_name;

 public:
  bool m_fill_clusters;
  bool m_fill_clusverbose;  // unmatched Svtx tracks
  bool m_fill_SvUnmatched;  // unmatched Svtx tracks
 private:
  // Tree Branch members:
  int nevent{-1};
  int nphg4{0};
  int nsvtx{0};
  int ntrackmatches{0};
  int nphg4_part{0};
  float centrality{0.};

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
  TH2 *h2_G4_nPixelsPhi;
  TH2 *h2_G4_nPixelsZ;
  TH2 *h2_Sv_nPixelsPhi;
  TH2 *h2_Sv_nPixelsZ;

  // Track tree
  int b_trackid{};
  bool b_is_g4track{};
  bool b_is_Svtrack{};
  bool b_is_matched{};

  float b_trkpt{};
  float b_trkphi{};
  float b_trketa{};

  int b_nclus{};
  int b_nclustpc{};
  int b_nclusmvtx{};
  int b_nclusintt{};

  float b_matchrat{};
  float b_matchrat_intt{};
  float b_matchrat_mvtx{};
  float b_matchrat_tpc{};

  std::vector<bool> b_clusmatch{};
  std::vector<float> b_clus_x{};
  std::vector<float> b_clus_y{};
  std::vector<float> b_clus_z{};
  std::vector<float> b_clus_r{};
  std::vector<int> b_clus_layer{};
  /* std::vector<int>   b_clus_nphibins   {}; */
  /* std::vector<int>   b_clus_ntbins     {}; */

  std::vector<float> b_clus_lphi{};      // l for local surface
  std::vector<float> b_clus_lphisize{};  // phi scaled by dphi
  std::vector<float> b_clus_lz{};
  std::vector<float> b_clus_lzsize{};

  // Maps of hits within clusters
  /* using BinData = std::pair<int,float>;  // index + energy for a given bin (MVTX, INTT, or TPC) */
  using VecVecInt = std::vector<std::vector<int>>;
  /* using VecVecFloat = std::vector<std::vector<float>>; */
  VecVecInt b_phibins;
  VecVecInt b_phibins_cut;
  VecVecInt b_zbins;
  VecVecInt b_zbins_cut;

  VecVecInt b_phibinsE;
  VecVecInt b_phibinsE_cut;
  VecVecInt b_zbinsE;
  VecVecInt b_zbinsE_cut;

  void fill_cluster_branches(TrkrClusLoc &, bool isPHG4 = true);
  /* void pushback_verb_bins(TrkrDefs::cluskey, bool isPHG4=true); */
};

#endif  // G4EVAL_FILLCLUSMATCHTREE_H

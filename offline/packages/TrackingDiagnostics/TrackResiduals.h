// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRESIDUALS_H
#define TRACKRESIDUALS_H

#include <tpc/TpcClusterMover.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

class TrkrCluster;
class PHCompositeNode;
class ActsGeometry;
class SvtxTrack;
class TrackSeed;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeomContainer;
class PHG4CylinderGeomContainer;
class TpcDistortionCorrectionContainer;
class TrackResiduals : public SubsysReco
{
 public:
  TrackResiduals(const std::string &name = "TrackResiduals");

  ~TrackResiduals() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void outfileName(const std::string &name) { m_outfileName = name; }
  void alignment(bool align) { m_doAlignment = align; }
  void alignmentmapName(const std::string &name) { m_alignmentMapName = name; }
  void trackmapName(const std::string &name) { m_trackMapName = name; }
  void clusterTree() { m_doClusters = true; }
  void vertexTree() { m_doVertex = true; }
  void hitTree() { m_doHits = true; }
  void eventTree() { m_doEventTree = true; }
  void MatchedTracksOnly() { m_doMatchedOnly = true; }
  void ppmode() { m_ppmode = true; }
  void convertSeeds(bool flag) { m_convertSeeds = flag; }
  void dropClustersNoState(bool flag) { m_dropClustersNoState = flag; }
  void zeroField() { m_zeroField = true; }
  void runnumber(const int run) { m_runnumber = run; }
  void segment(const int seg) { m_segment = seg; }
  void linefitAll() { m_linefitTPCOnly = false; }
  void setClusterMinSize(unsigned int size) { m_min_cluster_size = size; }
  void failedTree() { m_doFailedSeeds = true; }
  void setSegment(const int segment) { m_segment = segment; }

  void set_doMicromegasOnly(bool value) { m_doMicromegasOnly = value; }
  void setTrkrClusterContainerName(std::string &name) { m_clusterContainerName = name; }

  void set_use_clustermover(bool flag) { m_use_clustermover = flag; }

 private:
  void fillStatesWithLineFit(const TrkrDefs::cluskey &ckey,
                             TrkrCluster *cluster, ActsGeometry *geometry);
  void clearClusterStateVectors();
  void createBranches();
  static float convertTimeToZ(ActsGeometry *geometry, TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);
  void fillEventTree(PHCompositeNode *topNode);
  void fillClusterTree(TrkrClusterContainer *clusters, ActsGeometry *geometry);
  void fillHitTree(TrkrHitSetContainer *hitmap, ActsGeometry *geometry,
                   PHG4TpcCylinderGeomContainer *tpcGeom, PHG4CylinderGeomContainer *mvtxGeom,
                   PHG4CylinderGeomContainer *inttGeom, PHG4CylinderGeomContainer *mmGeom);
  void fillResidualTreeKF(PHCompositeNode *topNode);
  void fillResidualTreeSeeds(PHCompositeNode *topNode);
  void fillClusterBranchesKF(TrkrDefs::cluskey ckey, SvtxTrack *track,
                             const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> &global_moved,
                             PHCompositeNode *topNode);
  void fillClusterBranchesSeeds(TrkrDefs::cluskey ckey,  // SvtxTrack* track,
                                const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> &global,
                                PHCompositeNode *topNode);
  void lineFitClusters(std::vector<TrkrDefs::cluskey> &keys, TrkrClusterContainer *clusters, const short int &crossing);
  void circleFitClusters(std::vector<TrkrDefs::cluskey> &keys, TrkrClusterContainer *clusters, const short int &crossing);
  void fillStatesWithCircleFit(const TrkrDefs::cluskey &key, TrkrCluster *cluster,
                               Acts::Vector3 &glob, ActsGeometry *geometry);
  void fillVertexTree(PHCompositeNode *topNode);
  void fillFailedSeedTree(PHCompositeNode *topNode, std::set<unsigned int> &tpc_seed_ids);
  static float calc_dedx(TrackSeed *tpcseed, TrkrClusterContainer *clustermap, PHG4TpcCylinderGeomContainer *tpcGeom);

  bool m_use_clustermover = true;

  std::string m_outfileName = "";
  TFile *m_outfile = nullptr;
  TTree *m_tree = nullptr;
  TTree *m_clustree = nullptr;
  TTree *m_eventtree = nullptr;
  TTree *m_hittree = nullptr;
  TTree *m_vertextree = nullptr;
  TTree *m_failedfits = nullptr;

  bool m_doVertex = false;
  bool m_doClusters = false;
  bool m_doHits = false;
  bool m_doEventTree = false;
  bool m_zeroField = false;
  bool m_doFailedSeeds = false;
  bool m_doMatchedOnly = false;

  TpcClusterMover m_clusterMover;
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  ClusterErrorPara m_clusErrPara;
  std::string m_alignmentMapName = "SvtxAlignmentStateMap";
  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_clusterContainerName = "TRKR_CLUSTER";

  bool m_doAlignment = false;
  bool m_ppmode = false;
  bool m_convertSeeds = false;
  bool m_linefitTPCOnly = true;
  bool m_dropClustersNoState = false;
  unsigned int m_min_cluster_size = 0;

  bool m_doMicromegasOnly = false;

  int m_event = 0;
  int m_segment = std::numeric_limits<int>::quiet_NaN();
  int m_runnumber = std::numeric_limits<int>::quiet_NaN();
  int m_ntpcclus = std::numeric_limits<int>::quiet_NaN();
  float m_totalmbd = std::numeric_limits<float>::quiet_NaN();
  std::vector<int> m_firedTriggers;
  uint64_t m_gl1BunchCrossing = std::numeric_limits<uint64_t>::quiet_NaN();

  //! Event level quantities
  int m_nmvtx_all = std::numeric_limits<int>::quiet_NaN();
  int m_nintt_all = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc_hits0 = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc_hits1 = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc_clus0 = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc_clus1 = std::numeric_limits<int>::quiet_NaN();
  int m_nmms_all = std::numeric_limits<int>::quiet_NaN();
  int m_nsiseed = std::numeric_limits<int>::quiet_NaN();
  int m_ntpcseed = std::numeric_limits<int>::quiet_NaN();
  int m_ntracks_all = std::numeric_limits<int>::quiet_NaN();
  std::vector<int> m_ntpc_clus_sector;

  //! Track level quantities
  uint64_t m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
  unsigned int m_trackid = std::numeric_limits<unsigned int>::quiet_NaN();
  int m_crossing = std::numeric_limits<int>::quiet_NaN();
  int m_crossing_estimate = std::numeric_limits<int>::quiet_NaN();
  unsigned int m_tpcid = std::numeric_limits<unsigned int>::quiet_NaN();
  unsigned int m_silid = std::numeric_limits<unsigned int>::quiet_NaN();
  float m_px = std::numeric_limits<float>::quiet_NaN();
  float m_py = std::numeric_limits<float>::quiet_NaN();
  float m_pz = std::numeric_limits<float>::quiet_NaN();
  float m_pt = std::numeric_limits<float>::quiet_NaN();
  float m_eta = std::numeric_limits<float>::quiet_NaN();
  float m_phi = std::numeric_limits<float>::quiet_NaN();
  float m_deltapt = std::numeric_limits<float>::quiet_NaN();
  int m_charge = std::numeric_limits<int>::quiet_NaN();
  float m_quality = std::numeric_limits<float>::quiet_NaN();
  float m_chisq = std::numeric_limits<float>::quiet_NaN();
  float m_ndf = std::numeric_limits<float>::quiet_NaN();
  int m_nhits = std::numeric_limits<int>::quiet_NaN();
  int m_nmaps = std::numeric_limits<int>::quiet_NaN();
  int m_nmapsstate = std::numeric_limits<int>::quiet_NaN();
  int m_nintt = std::numeric_limits<int>::quiet_NaN();
  int m_ninttstate = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc = std::numeric_limits<int>::quiet_NaN();
  int m_ntpcstate = std::numeric_limits<int>::quiet_NaN();
  int m_nmms = std::numeric_limits<int>::quiet_NaN();
  int m_nmmsstate = std::numeric_limits<int>::quiet_NaN();
  unsigned int m_vertexid = std::numeric_limits<unsigned int>::quiet_NaN();
  int m_vertex_crossing = std::numeric_limits<int>::quiet_NaN();
  float m_vx = std::numeric_limits<float>::quiet_NaN();
  float m_vy = std::numeric_limits<float>::quiet_NaN();
  float m_vz = std::numeric_limits<float>::quiet_NaN();
  int m_vertex_ntracks = std::numeric_limits<int>::quiet_NaN();
  float m_pcax = std::numeric_limits<float>::quiet_NaN();
  float m_pcay = std::numeric_limits<float>::quiet_NaN();
  float m_pcaz = std::numeric_limits<float>::quiet_NaN();
  float m_rzslope = std::numeric_limits<float>::quiet_NaN();
  float m_rzint = std::numeric_limits<float>::quiet_NaN();
  float m_xyslope = std::numeric_limits<float>::quiet_NaN();
  float m_xyint = std::numeric_limits<float>::quiet_NaN();
  float m_yzslope = std::numeric_limits<float>::quiet_NaN();
  float m_yzint = std::numeric_limits<float>::quiet_NaN();
  float m_R = std::numeric_limits<float>::quiet_NaN();
  float m_X0 = std::numeric_limits<float>::quiet_NaN();
  float m_Y0 = std::numeric_limits<float>::quiet_NaN();
  float m_dcaxy = std::numeric_limits<float>::quiet_NaN();
  float m_dcaz = std::numeric_limits<float>::quiet_NaN();
  float m_tracklength = std::numeric_limits<float>::quiet_NaN();

  float m_silseedx = std::numeric_limits<float>::quiet_NaN();
  float m_silseedy = std::numeric_limits<float>::quiet_NaN();
  float m_silseedz = std::numeric_limits<float>::quiet_NaN();
  float m_silseedpx = std::numeric_limits<float>::quiet_NaN();
  float m_silseedpy = std::numeric_limits<float>::quiet_NaN();
  float m_silseedpz = std::numeric_limits<float>::quiet_NaN();
  int m_silseedcharge = std::numeric_limits<int>::quiet_NaN();
  float m_silseedphi = std::numeric_limits<float>::quiet_NaN();
  float m_silseedeta = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedx = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedy = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedz = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedpx = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedpy = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedpz = std::numeric_limits<float>::quiet_NaN();
  int m_tpcseedcharge = std::numeric_limits<int>::quiet_NaN();
  float m_tpcseedphi = std::numeric_limits<float>::quiet_NaN();
  float m_tpcseedeta = std::numeric_limits<float>::quiet_NaN();

  float m_dedx = std::numeric_limits<float>::quiet_NaN();

  //! hit tree info
  uint32_t m_hitsetkey = std::numeric_limits<uint32_t>::quiet_NaN();
  float m_hitgx = std::numeric_limits<float>::quiet_NaN();
  float m_hitgy = std::numeric_limits<float>::quiet_NaN();
  float m_hitgz = std::numeric_limits<float>::quiet_NaN();
  int m_hitlayer = std::numeric_limits<int>::quiet_NaN();
  int m_sector = std::numeric_limits<int>::quiet_NaN();
  int m_hitpad = std::numeric_limits<int>::quiet_NaN();
  int m_hittbin = std::numeric_limits<int>::quiet_NaN();
  int m_col = std::numeric_limits<int>::quiet_NaN();
  int m_row = std::numeric_limits<int>::quiet_NaN();
  int m_strip = std::numeric_limits<int>::quiet_NaN();
  float m_zdriftlength = std::numeric_limits<float>::quiet_NaN();
  
  float m_mbdvtxz = std::numeric_limits<float>::quiet_NaN();

  int m_ntracks = std::numeric_limits<int>::quiet_NaN();
  int m_nvertices = std::numeric_limits<int>::quiet_NaN();

  //! cluster tree info
  float m_sclusgr = std::numeric_limits<float>::quiet_NaN();
  float m_sclusphi = std::numeric_limits<float>::quiet_NaN();
  float m_scluseta = std::numeric_limits<float>::quiet_NaN();
  float m_adc = std::numeric_limits<float>::quiet_NaN();
  float m_clusmaxadc = std::numeric_limits<float>::quiet_NaN();
  int m_phisize = std::numeric_limits<int>::quiet_NaN();
  int m_zsize = std::numeric_limits<int>::quiet_NaN();
  float m_scluslx = std::numeric_limits<float>::quiet_NaN();
  float m_scluslz = std::numeric_limits<float>::quiet_NaN();
  float m_sclusgx = std::numeric_limits<float>::quiet_NaN();
  float m_sclusgy = std::numeric_limits<float>::quiet_NaN();
  float m_sclusgz = std::numeric_limits<float>::quiet_NaN();
  int m_scluslayer = std::numeric_limits<int>::quiet_NaN();
  float m_scluselx = std::numeric_limits<float>::quiet_NaN();
  float m_scluselz = std::numeric_limits<float>::quiet_NaN();
  int m_clussector = std::numeric_limits<int>::quiet_NaN();
  int m_side = std::numeric_limits<int>::quiet_NaN();
  int m_staveid = std::numeric_limits<int>::quiet_NaN();
  int m_chipid = std::numeric_limits<int>::quiet_NaN();
  int m_strobeid = std::numeric_limits<int>::quiet_NaN();
  int m_ladderzid = std::numeric_limits<int>::quiet_NaN();
  int m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
  int m_timebucket = std::numeric_limits<int>::quiet_NaN();
  int m_segtype = std::numeric_limits<int>::quiet_NaN();
  int m_tileid = std::numeric_limits<int>::quiet_NaN();

  //! clusters on track information
  std::vector<float> m_clusAdc;
  std::vector<float> m_clusMaxAdc;
  std::vector<float> m_cluslx;
  std::vector<float> m_cluslz;
  std::vector<float> m_cluselx;
  std::vector<float> m_cluselz;
  std::vector<float> m_clusgx;
  std::vector<float> m_clusgy;
  std::vector<float> m_clusgz;
  std::vector<float> m_clusgxunmoved;
  std::vector<float> m_clusgyunmoved;
  std::vector<float> m_clusgzunmoved;
  std::vector<float> m_clusgr;
  std::vector<int> m_clsector;
  std::vector<int> m_clside;
  std::vector<int> m_cluslayer;
  std::vector<int> m_clusphisize;
  std::vector<int> m_cluszsize;
  std::vector<int> m_clusedge;
  std::vector<int> m_clusoverlap;
  std::vector<uint64_t> m_cluskeys;
  std::vector<float> m_idealsurfcenterx;
  std::vector<float> m_idealsurfcentery;
  std::vector<float> m_idealsurfcenterz;
  std::vector<float> m_idealsurfnormx;
  std::vector<float> m_idealsurfnormy;
  std::vector<float> m_idealsurfnormz;
  std::vector<float> m_missurfcenterx;
  std::vector<float> m_missurfcentery;
  std::vector<float> m_missurfcenterz;
  std::vector<float> m_missurfnormx;
  std::vector<float> m_missurfnormy;
  std::vector<float> m_missurfnormz;
  std::vector<float> m_clusgxideal;
  std::vector<float> m_clusgyideal;
  std::vector<float> m_clusgzideal;
  std::vector<float> m_idealsurfalpha;
  std::vector<float> m_idealsurfbeta;
  std::vector<float> m_idealsurfgamma;
  std::vector<float> m_missurfalpha;
  std::vector<float> m_missurfbeta;
  std::vector<float> m_missurfgamma;

  //! states on track information
  std::vector<float> m_statelx;
  std::vector<float> m_statelz;
  std::vector<float> m_stateelx;
  std::vector<float> m_stateelz;
  std::vector<float> m_stategx;
  std::vector<float> m_stategy;
  std::vector<float> m_stategz;
  std::vector<float> m_statepx;
  std::vector<float> m_statepy;
  std::vector<float> m_statepz;
  std::vector<float> m_statepl;

  std::vector<float> m_statelxglobderivdx;
  std::vector<float> m_statelxglobderivdy;
  std::vector<float> m_statelxglobderivdz;
  std::vector<float> m_statelxglobderivdalpha;
  std::vector<float> m_statelxglobderivdbeta;
  std::vector<float> m_statelxglobderivdgamma;

  std::vector<float> m_statelxlocderivd0;
  std::vector<float> m_statelxlocderivz0;
  std::vector<float> m_statelxlocderivphi;
  std::vector<float> m_statelxlocderivtheta;
  std::vector<float> m_statelxlocderivqop;

  std::vector<float> m_statelzglobderivdx;
  std::vector<float> m_statelzglobderivdy;
  std::vector<float> m_statelzglobderivdz;
  std::vector<float> m_statelzglobderivdalpha;
  std::vector<float> m_statelzglobderivdbeta;
  std::vector<float> m_statelzglobderivdgamma;

  std::vector<float> m_statelzlocderivd0;
  std::vector<float> m_statelzlocderivz0;
  std::vector<float> m_statelzlocderivphi;
  std::vector<float> m_statelzlocderivtheta;
  std::vector<float> m_statelzlocderivqop;
};

#endif  // TRACKRESIDUALS_H

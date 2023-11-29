// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRESIDUALS_H
#define TRACKRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <string>

#include <trackbase/TrkrDefs.h>
#include <trackbase/ClusterErrorPara.h>

#include <cmath>
#include <iostream>
#include <limits>

class TrkrCluster;
class PHCompositeNode;
class ActsGeometry;
class SvtxTrack;

class TrackResiduals : public SubsysReco
{
 public:
  TrackResiduals(const std::string &name = "TrackResiduals");

  ~TrackResiduals() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void outfileName(const std::string &name) { m_outfileName = name; }
  void alignment(bool align) { m_doAlignment = align; }
  void alignmentmapName(const std::string &name) { m_alignmentMapName = name; }
  void trackmapName(const std::string &name) { m_trackMapName = name; }

 private:
  void clearClusterStateVectors();
  void createBranches();
  float convertTimeToZ(ActsGeometry *geometry, TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);

  void fillClusterBranches(TrkrDefs::cluskey ckey, SvtxTrack *track,
                           PHCompositeNode *topNode);

  std::string m_outfileName = "";
  TFile *m_outfile = nullptr;
  TTree *m_tree = nullptr;

  ClusterErrorPara m_clusErrPara;
  std::string m_alignmentMapName = "SvtxAlignmentStateMap";
  std::string m_trackMapName = "SvtxTrackMap";

  bool m_doAlignment = false;

  int m_event = 0;
  //! Track level quantities
  unsigned int m_trackid = std::numeric_limits<unsigned int>::quiet_NaN();
  unsigned int m_crossing = std::numeric_limits<unsigned int>::quiet_NaN();
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
  int m_nintt = std::numeric_limits<int>::quiet_NaN();
  int m_ntpc = std::numeric_limits<int>::quiet_NaN();
  int m_nmms = std::numeric_limits<int>::quiet_NaN();
  unsigned int m_vertexid = std::numeric_limits<unsigned int>::quiet_NaN();
  float m_vx = std::numeric_limits<float>::quiet_NaN();
  float m_vy = std::numeric_limits<float>::quiet_NaN();
  float m_vz = std::numeric_limits<float>::quiet_NaN();
  float m_pcax = std::numeric_limits<float>::quiet_NaN();
  float m_pcay = std::numeric_limits<float>::quiet_NaN();
  float m_pcaz = std::numeric_limits<float>::quiet_NaN();

  //! clusters on track information
  std::vector<float> m_cluslx;
  std::vector<float> m_cluslz;
  std::vector<float> m_cluselx;
  std::vector<float> m_cluselz;
  std::vector<float> m_clusgx;
  std::vector<float> m_clusgy;
  std::vector<float> m_clusgz;
  std::vector<int> m_cluslayer;
  std::vector<int> m_clussize;
  std::vector<int> m_clusedge;
  std::vector<int> m_clusoverlap;
  std::vector<uint32_t> m_clushitsetkey;
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

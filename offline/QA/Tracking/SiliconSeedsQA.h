// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SILICONSEEDSQA_H
#define SILICONSEEDSQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class SvtxTrack;
class PHCompositeNode;
class TH1;
class TH2;
class TProfile2D;

class SiliconSeedsQA : public SubsysReco
{
 public:
  SiliconSeedsQA(const std::string &name = "SiliconSeedsQA");

  ~SiliconSeedsQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;
  void setTrackMapName(const std::string &name) { m_trackMapName = name; }
  void setVertexMapName(const std::string &name) { m_vertexMapName = name; }

 private:
  static std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack *track);
  void createHistos();
  std::string getHistoPrefix() const;

  std::string m_clusterContainerName = "TRKR_CLUSTER";
  std::string m_actsgeometryName = "ActsGeometry";
  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";

  TH1 *h_ntrack1d = nullptr;
  TH2 *h_ntrack = nullptr;
  TH1 *h_nmaps = nullptr;
  TH1 *h_nintt = nullptr;
  TH2 *h_nmaps_nintt = nullptr;
  TProfile2D *h_avgnclus_eta_phi = nullptr;
  TH1 *h_trackcrossing = nullptr;
  TH1 *h_trackchi2ndf = nullptr;
  TH2 *h_dcaxyorigin_phi = nullptr;
  TH2 *h_dcaxyvtx_phi = nullptr;
  TH2 *h_dcazorigin_phi = nullptr;
  TH2 *h_dcazvtx_phi = nullptr;
  TH1 *h_ntrack_isfromvtx = nullptr;
  TH1 *h_trackpt_inclusive = nullptr;
  TH1 *h_trackpt_pos = nullptr;
  TH1 *h_trackpt_neg = nullptr;
  TH1 *h_ntrack_IsPosCharge = nullptr;

  TH1 *h_nvertex = nullptr;
  TH1 *h_vx = nullptr;
  TH1 *h_vy = nullptr;
  TH1 *h_vz = nullptr;
  TH2 *h_vx_vy = nullptr;
  TH1 *h_vt = nullptr;
  TH1 *h_vcrossing = nullptr;
  TH1 *h_vchi2dof = nullptr;
  TH1 *h_ntrackpervertex = nullptr;
};

#endif  // SILICONSEEDSQA_H

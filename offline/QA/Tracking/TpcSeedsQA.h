// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCSEEDSQA_H
#define TPCSEEDSQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <globalvertex/SvtxVertexMap.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <string>
#include <vector>

class SvtxTrack;
class PHCompositeNode;

class TpcSeedsQA : public SubsysReco
{
 public:
  TpcSeedsQA(const std::string &name = "TpcSeedsQA");

  ~TpcSeedsQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

 private:
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack *track);
  void createHistos();
  std::string getHistoPrefix() const;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";

  TrkrClusterContainer* clustermap;
  ActsGeometry* geometry;
  SvtxTrackMap* trackmap;
  SvtxVertexMap* vertexmap;

  TH1* h_ntrack1d;
  TH1* h_ntrack1d_pos;
  TH1* h_ntrack1d_neg;
  TH2* h_ntrack_pos;
  TH2* h_ntrack_neg;
  TH1* h_ntpc_pos;
  TH1* h_ntpc_neg;
  TH1* h_ntpot_pos;
  TH1* h_ntpot_neg;
  TProfile2D* h_avgnclus_eta_phi_pos;
  TProfile2D* h_avgnclus_eta_phi_neg;
  TH1* h_trackcrossing_pos;
  TH1* h_trackcrossing_neg;
  TH2* h_dcaxyorigin_phi_pos;
  TH2* h_dcaxyorigin_phi_neg;
  TH2* h_dcaxyvtx_phi_pos;
  TH2* h_dcaxyvtx_phi_neg;
  TH2* h_dcazorigin_phi_pos;
  TH2* h_dcazorigin_phi_neg;
  TH2* h_dcazvtx_phi_pos;
  TH2* h_dcazvtx_phi_neg;
  TH1* h_ntrack_isfromvtx_pos;
  TH1* h_ntrack_isfromvtx_neg;
  TH1* h_cluster_phisize1_fraction_pos;
  TH1* h_cluster_phisize1_fraction_neg;

  TH1* h_nvertex;
  TH1* h_vx;
  TH1* h_vy;
  TH2* h_vx_vy;
  TH1* h_vz;
  TH1* h_vt;
  TH1* h_vcrossing;
  TH1* h_vchi2dof;
  TH1* h_ntrackpervertex;

};

#endif  // TPCSEEDSQA_H

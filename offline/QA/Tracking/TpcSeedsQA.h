// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCSEEDSQA_H
#define TPCSEEDSQA_H

#include <fun4all/SubsysReco.h>
/* #include <trackbase/ActsGeometry.h> */
#include <trackbase/TrkrDefs.h>

/* #include <trackbase/TrkrClusterContainer.h> */
/* #include <trackbase/TrkrCluster.h> */
/* #include <trackbase_historic/SvtxTrackMap.h> */
/* #include <globalvertex/SvtxVertexMap.h> */
/* #include <TH1.h> */
/* #include <TH2.h> */
/* #include <TProfile2D.h> */

#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TH1;
class TH2;
class TrkrClusterContainer;
class TProfile2D;
class SvtxVertexMap;

class TpcSeedsQA : public SubsysReco
{
 public:
  TpcSeedsQA(const std::string& name = "TpcSeedsQA");

  ~TpcSeedsQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;

  int End(PHCompositeNode* topNode) override;
  void setTrackMapName(const std::string& name) { m_trackMapName = name; }
  void setVertexMapName(const std::string& name) { m_vertexMapName = name; }

 private:
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  void createHistos();
  std::string getHistoPrefix() const;

  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_vertexMapName{"SvtxVertexMap"};

  TrkrClusterContainer* clustermap{nullptr};
  ActsGeometry* geometry{nullptr};
  SvtxTrackMap* trackmap{nullptr};
  SvtxVertexMap* vertexmap{nullptr};

  TH1* h_ntrack1d{nullptr};
  TH1* h_ntrack1d_pos{nullptr};
  TH1* h_ntrack1d_neg{nullptr};
  TH2* h_ntrack_pos{nullptr};
  TH2* h_ntrack_neg{nullptr};
  TH1* h_ntpc_pos{nullptr};
  TH1* h_ntpc_neg{nullptr};
  TH1* h_ntpot_pos{nullptr};
  TH1* h_ntpot_neg{nullptr};
  TProfile2D* h_avgnclus_eta_phi_pos{nullptr};
  TProfile2D* h_avgnclus_eta_phi_neg{nullptr};
  TH1* h_trackcrossing_pos{nullptr};
  TH1* h_trackcrossing_neg{nullptr};
  TH2* h_dcaxyorigin_phi_pos{nullptr};
  TH2* h_dcaxyorigin_phi_neg{nullptr};
  TH2* h_dcaxyvtx_phi_pos{nullptr};
  TH2* h_dcaxyvtx_phi_neg{nullptr};
  TH2* h_dcazorigin_phi_pos{nullptr};
  TH2* h_dcazorigin_phi_neg{nullptr};
  TH2* h_dcazvtx_phi_pos{nullptr};
  TH2* h_dcazvtx_phi_neg{nullptr};
  TH1* h_ntrack_isfromvtx_pos{nullptr};
  TH1* h_ntrack_isfromvtx_neg{nullptr};
  TH1* h_cluster_phisize1_fraction_pos{nullptr};
  TH1* h_cluster_phisize1_fraction_neg{nullptr};

  TH1* h_nvertex{nullptr};
  TH1* h_vx{nullptr};
  TH1* h_vy{nullptr};
  TH2* h_vx_vy{nullptr};
  TH1* h_vz{nullptr};
  TH1* h_vt{nullptr};
  TH1* h_vcrossing{nullptr};
  TH1* h_vchi2dof{nullptr};
  TH1* h_ntrackpervertex{nullptr};
};

#endif  // TPCSEEDSQA_H

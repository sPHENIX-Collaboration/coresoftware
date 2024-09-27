// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCSEEDSQA_H
#define TPCSEEDSQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TH1;
class TH2;
class TrkrClusterContainer;
class TProfile;
class TProfile2D;
class TNtuple;
class SvtxVertexMap;
class TrackSeedContainer;
class PHG4TpcCylinderGeomContainer;
class TrackSeed;

class TpcSeedsQA : public SubsysReco
{
 public:
  TpcSeedsQA(const std::string& name = "TpcSeedsQA");

  ~TpcSeedsQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode* topNode) override;

  void setClusterContainerName(const std::string& name) { m_clusterContainerName = name; }
  std::string getClusterContainerName() { return m_clusterContainerName; }
  void setActsGeomName(const std::string& name) { m_actsGeomName = name; }
  std::string getActsGeomName() { return m_actsGeomName; }
  void setG4GeomName(const std::string& name) { m_g4GeomName = name; }
  std::string getG4GeomName() { return m_g4GeomName; }
  void setTrackMapName(const std::string& name) { m_trackMapName = name; }
  std::string getTrackMapName() { return m_trackMapName; }
  void setVertexMapName(const std::string& name) { m_vertexMapName = name; }
  std::string gettVertexMapName() { return m_vertexMapName; }
  float calc_dedx(TrackSeed* tpcseed);
  void setSegment(const int segment) { m_segment = segment; }
  void segment(const int seg) { m_segment = seg; }

 private:
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  void createHistos();
  std::string getHistoPrefix() const;
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;
  std::pair<float, float> cal_tpc_eta_min_max(float vtxz);

  std::string m_clusterContainerName{"TRKR_CLUSTER"};
  std::string m_actsGeomName{"ActsGeometry"};
  std::string m_g4GeomName{"CYLINDERCELLGEOM_SVTX"};
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_vertexMapName{"SvtxVertexMap"};

  TrkrClusterContainer* clustermap{nullptr};
  ActsGeometry* actsgeom{nullptr};
  PHG4TpcCylinderGeomContainer* g4geom{nullptr};
  SvtxTrackMap* trackmap{nullptr};
  SvtxVertexMap* vertexmap{nullptr};

  TH1* h_ntrack1d{nullptr};
  TH1* h_ntrack1d_pos{nullptr};
  TH1* h_ntrack1d_neg{nullptr};
  TH1* h_ntrack1d_ptg1{nullptr};
  TH1* h_ntrack1d_ptg1_pos{nullptr};
  TH1* h_ntrack1d_ptg1_neg{nullptr};
  TH1* h_pt{nullptr};
  TH1* h_pt_pos{nullptr};
  TH1* h_pt_neg{nullptr};
  TH2* h_ntrack_pos{nullptr};
  TH2* h_ntrack_neg{nullptr};
  TH1* h_ntpc_fullpt_pos{nullptr};
  TH1* h_ntpc_fullpt_neg{nullptr};
  TH1* h_ntpc_pos{nullptr};
  TH1* h_ntpc_neg{nullptr};
  TH2* h_ntpc_quality_pos{nullptr};
  TH2* h_ntpc_quality_neg{nullptr};
  TH1* h_ntpot_pos{nullptr};
  TH1* h_ntpot_neg{nullptr};
  TProfile2D* h_avgnclus_eta_phi_pos{nullptr};
  TProfile2D* h_avgnclus_eta_phi_neg{nullptr};
  // TH1* h_trackcrossing_pos{nullptr};
  // TH1* h_trackcrossing_neg{nullptr};
  TH2* h_dcaxyorigin_phi_north_pos{nullptr};
  TH2* h_dcaxyorigin_phi_south_pos{nullptr};
  TH2* h_dcaxyorigin_phi_north_neg{nullptr};
  TH2* h_dcaxyorigin_phi_south_neg{nullptr};
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
  // TH1* h_vcrossing{nullptr};
  TH1* h_vchi2dof{nullptr};
  TH1* h_ntrackpervertex{nullptr};

  TH2* h_dedx{nullptr};
  TH1* h_mip_dedx{nullptr};

  TH2* h_adc_sector[3] = {nullptr};
  TProfile* h_onepad_frac[3] = {nullptr};

  TH1* h_cluster_phisize1_fraction_side0[3] = {nullptr};
  TH1* h_cluster_phisize1_fraction_side1[3] = {nullptr};

  TH1* h_clusphisize1pt_side0[3] = {nullptr};
  TH1* h_clusphisize1pt_side1[3] = {nullptr};
  TH1* h_clusphisizegeq1pt_side0[3] = {nullptr};
  TH1* h_clusphisizegeq1pt_side1[3] = {nullptr};

  TH2* h_cluster_phisize1_fraction_pt_side0[3] = {nullptr};
  TH2* h_cluster_phisize1_fraction_pt_side1[3] = {nullptr};

  TH1* h_cluster_phisize1_fraction_mean_side0[3] = {nullptr};
  TH1* h_cluster_phisize1_fraction_mean_side1[3] = {nullptr};

  TH1* h_cluster_phisize1_fraction_mean_numerator_side0[3] = {nullptr};
  TH1* h_cluster_phisize1_fraction_mean_numerator_side1[3] = {nullptr};

  TH1* h_cluster_phisize1_fraction_mean_denominator_side0[3] = {nullptr};
  TH1* h_cluster_phisize1_fraction_mean_denominator_side1[3] = {nullptr};

  TNtuple* nt_sector_event_summary = {nullptr};

  int m_bco = 0;;
  int m_event = 0;
  int m_segment = 0;

  double frac_side0_pt[3][4] = {{0}};
  double frac_side1_pt[3][4] = {{0}};

  double num_track_side0_pt[3][4] = {{0}};
  double num_track_side1_pt[3][4] = {{0}};

  //! global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  float m_px = std::numeric_limits<float>::quiet_NaN();
  float m_py = std::numeric_limits<float>::quiet_NaN();
  float m_pz = std::numeric_limits<float>::quiet_NaN();
  float m_pt = std::numeric_limits<float>::quiet_NaN();
  float m_ptot = std::numeric_limits<float>::quiet_NaN();
  float m_charge = std::numeric_limits<float>::quiet_NaN();
  float m_dedx = std::numeric_limits<float>::quiet_NaN();

  int m_ntpc = std::numeric_limits<int>::quiet_NaN();
  std::vector<float> m_clusgz;
  std::vector<int> m_cluslayer;
  std::vector<int> m_clusphisize;
  std::vector<int> m_cluszsize;
  std::vector<int> m_region;

  struct PhiHistoList
  {
    TH1* cphisize1pT_side0 = nullptr;
    TH1* cphisize1pT_side1 = nullptr;
    TH1* cphisizegeq1pT_side0 = nullptr;
    TH1* cphisizegeq1pT_side1 = nullptr;
    int ntpc_side0 = 0;
    int ntpc_side0_phisize1 = 0;
    int ntpc_side1 = 0;
    int ntpc_side1_phisize1 = 0;

    void Clear()
    {
      ntpc_side0 = 0;
      ntpc_side0_phisize1 = 0;
      ntpc_side1 = 0;
      ntpc_side1_phisize1 = 0;
    }
  };

  using PhiHistoMap = std::map<int, PhiHistoList>;
  PhiHistoMap phihistos;
};

#endif  // TPCSEEDSQA_H

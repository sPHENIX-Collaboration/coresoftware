// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef STATECLUSTERRESIDUALSQA_H
#define STATECLUSTERRESIDUALSQA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <limits>
#include <cmath>
#include <cfloat>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;

struct ResidualHistConfig
{
  std::string name  = "h_StateClusterResidualsQA_";
  std::string title = ";Residual [cm];Entries";

  int min_mvtx_clusters = 0;
  int max_mvtx_clusters = 3;
  int min_intt_clusters = 0;
  int max_intt_clusters = 4;
  int min_tpc_clusters = 0;
  int max_tpc_clusters = 48;

  float phi_min = -M_PI;
  float phi_max =  M_PI;
  float eta_min = -1.1;
  float eta_max = 1.1;

  float pt_min = 0.0;
  float pt_max = FLT_MAX;

  int charge = 0;

  float rphi_local_lower = -0.5;
  float rphi_local_upper = 0.5;
  float z_local_lower = -0.5;
  float z_local_upper = 0.5;
  float x_lower = -0.5;
  float x_upper = 0.5;
  float y_lower = -0.5;
  float y_upper = 0.5;
  float z_lower = -0.5;
  float z_upper = 0.5;
};

class StateClusterResidualsQA : public SubsysReco
{
 public:
  StateClusterResidualsQA(const std::string& name = "StateClusterResidualsQA");
  ~StateClusterResidualsQA() override = default;

  /// sets the name of node to retrieve the track map from (default member value is "SvtxTrackMap")
  void set_track_map_name(std::string const& track_map_node_name) { m_track_map_node_name = track_map_node_name; }

  StateClusterResidualsQA& addHistogram(const std::string& name)
  {
    ResidualHistConfig cfg;
    cfg.name += name;
    m_pending.push_back(cfg);
    return *this;
  }
  StateClusterResidualsQA& setNMvtx(int min, int max)
  {
    m_pending.back().min_mvtx_clusters = min;
    m_pending.back().max_mvtx_clusters = max;
    return *this;
  }
  StateClusterResidualsQA& setNIntt(int min, int max)
  {
    m_pending.back().min_intt_clusters = min;
    m_pending.back().max_intt_clusters = max;
    return *this;
  }
  StateClusterResidualsQA& setNTpc(int min, int max)
  {
    m_pending.back().min_tpc_clusters = min;
    m_pending.back().max_tpc_clusters = max;
    return *this;
  }
  StateClusterResidualsQA& setPhiRange(float min, float max)
  {
    m_pending.back().phi_min = min;
    m_pending.back().phi_max = max;
    return *this;
  }
  StateClusterResidualsQA& setEtaRange(float min, float max)
  {
    m_pending.back().eta_min = min;
    m_pending.back().eta_max = max;
    return *this;
  }
  StateClusterResidualsQA& setPtRange(float min, float max)
  {
    m_pending.back().pt_min = min;
    m_pending.back().pt_max = max;
    return *this;
  }
  StateClusterResidualsQA& setXRange(float min, float max)
  {
    m_pending.back().x_lower = min;
    m_pending.back().x_upper = max;
    return *this;
  }
  StateClusterResidualsQA& setYRange(float min, float max)
  {
    m_pending.back().y_lower = min;
    m_pending.back().y_upper = max;
    return *this;
  }
  StateClusterResidualsQA& setZRange(float min, float max)
  {
    m_pending.back().z_lower = min;
    m_pending.back().z_upper = max;
    return *this;
  }
  StateClusterResidualsQA& setLocalRphiRange(float min, float max)
  {
    m_pending.back().rphi_local_lower = min;
    m_pending.back().rphi_local_upper = max;
    return *this;
  }
  StateClusterResidualsQA& setLocalZRange(float min, float max)
  {
    m_pending.back().z_local_lower = min;
    m_pending.back().z_local_upper = max;
    return *this;
  }
  StateClusterResidualsQA& setPositiveTracks()
  {
    m_pending.back().charge = 1;
    return *this;
  }
  StateClusterResidualsQA& setNegativeTracks()
  {
    m_pending.back().charge = -1;
    return *this;
  }
  
  void setUseLocalCoords()
  {
    m_use_local_coords = true;
  }

  void createHistos();

  int InitRun(PHCompositeNode*) override;

  int process_event(PHCompositeNode*) override;

  int EndRun(const int runnumber) override;
 
 private:
  std::vector<ResidualHistConfig> m_pending;

  std::string m_track_map_node_name = "SvtxTrackMap";
  std::string m_clusterContainerName = "TRKR_CLUSTER";

  int m_nBins = 50;
  bool m_use_local_coords = false;

  std::vector<TH1*> m_histograms_x{};
  std::vector<TH1*> m_histograms_y{};
  std::vector<TH1*> m_histograms_z{};
  std::vector<TH2*> m_histograms_layer_x{};
  std::vector<TH2*> m_histograms_layer_y{};
  std::vector<TH2*> m_histograms_layer_z{};
  std::vector<TH1*> m_histograms_phi_x{};
  std::vector<TH1*> m_histograms_phi_y{};
  std::vector<TH1*> m_histograms_phi_z{};
  std::vector<TH1*> m_histograms_eta_x{};
  std::vector<TH1*> m_histograms_eta_y{};
  std::vector<TH1*> m_histograms_eta_z{};
};

#endif  // TRACKFITTINGQA_H

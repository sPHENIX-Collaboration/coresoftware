// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONDISTORTIONS_H
#define QAG4SIMULATIONDISTORTIONS_H

#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <limits>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class TH1;
class TH2;
class TrkrClusterContainer;
class TTree;

class QAG4SimulationDistortions : public SubsysReco
{
 public:
  QAG4SimulationDistortions(const std::string &name = "QAG4SimulationDistortions");

  ~QAG4SimulationDistortions() override = default;

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

  //! track map name
  void set_trackmap_name(const std::string &value)
  {
    m_trackmapname = value;
  }

  void disableModuleEdgeCorr() { m_disable_module_edge_corr = true; }
  void disableStaticCorr() { m_disable_static_corr = true; }
  void disableAverageCorr() { m_disable_average_corr = true; }
  void disableFluctuationCorr() { m_disable_fluctuation_corr = true; }

  /// require micromegas to be present when extrapolating tracks to the TPC
  void setUseMicromegas(bool value)
  {
    m_useMicromegas = value;
  }

  /// set min track pt
  void setMinpT(double value)
  {
    m_minPT = value;
  }

 private:
  //! track map name
  std::string m_trackmapname = "SvtxSiliconMMTrackMap";

  std::string get_histo_prefix()
  {
    return std::string("h_") + Name() + std::string("_");
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack *track);
  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack *track);
  bool checkTrack(SvtxTrack *track);
  bool checkTPOTResidual(SvtxTrack *track);

  TH2 *h_beta{nullptr};
  TH2 *h_alpha{nullptr};
  TH2 *h_rphiResid{nullptr};
  TH2 *h_zResid{nullptr};
  TH2 *h_etaResid{nullptr};
  TH2 *h_etaResidLayer{nullptr};
  TH2 *h_zResidLayer{nullptr};
  TH2 *h_deltarphi_layer{nullptr};
  TH2 *h_deltaz_layer{nullptr};
  TH2 *h_statez_pulls{nullptr};
  TH2 *h_staterphi_pulls{nullptr};
  TH2 *h_clusz_pulls{nullptr};
  TH2 *h_clusrphi_pulls{nullptr};
  TH2 *h_nstates_vs_nclus{nullptr};

  TH1 *h_tpot_deltarphi{nullptr};
  TH1 *h_tpot_deltaz{nullptr};

  TTree *t_tree{nullptr};

  SvtxTrackMap *m_trackMap{nullptr};
  TrkrClusterContainer *m_clusterContainer{nullptr};
  ActsGeometry *m_tGeometry{nullptr};
  PHG4TpcGeomContainer *m_tpcGeom{nullptr};

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  int m_event{0};
  float m_tanAlpha{std::numeric_limits<float>::quiet_NaN()};
  float m_tanBeta{std::numeric_limits<float>::quiet_NaN()};
  float m_drphi{std::numeric_limits<float>::quiet_NaN()};
  float m_dz{std::numeric_limits<float>::quiet_NaN()};
  float m_clusR{std::numeric_limits<float>::quiet_NaN()};
  float m_clusPhi{std::numeric_limits<float>::quiet_NaN()};
  float m_clusZ{std::numeric_limits<float>::quiet_NaN()};
  float m_clusEta{std::numeric_limits<float>::quiet_NaN()};
  float m_stateR{std::numeric_limits<float>::quiet_NaN()};
  float m_statePhi{std::numeric_limits<float>::quiet_NaN()};
  float m_stateZ{std::numeric_limits<float>::quiet_NaN()};
  float m_stateEta{std::numeric_limits<float>::quiet_NaN()};
  float m_stateRPhiErr{std::numeric_limits<float>::quiet_NaN()};
  float m_stateZErr{std::numeric_limits<float>::quiet_NaN()};
  float m_clusRPhiErr{std::numeric_limits<float>::quiet_NaN()};
  float m_clusZErr{std::numeric_limits<float>::quiet_NaN()};

  int m_layer{-1};
  float m_statePt{std::numeric_limits<float>::quiet_NaN()};
  float m_statePz{std::numeric_limits<float>::quiet_NaN()};
  float m_trackPt{std::numeric_limits<float>::quiet_NaN()};
  float m_trackdEdx{std::numeric_limits<float>::quiet_NaN()};
  int m_charge{-10};
  int m_crossing{-10};

  TrkrDefs::cluskey m_cluskey{TrkrDefs::CLUSKEYMAX};

  /// disable distortion correction
  bool m_disable_module_edge_corr{false};
  bool m_disable_static_corr{false};
  bool m_disable_average_corr{false};
  bool m_disable_fluctuation_corr{false};

  ///@name counters
  //@{
  int m_total_tracks{0};
  int m_accepted_tracks{0};

  int m_total_states{0};
  int m_accepted_states{0};
  //@}

  bool m_useMicromegas{true};
  double m_minPT{0.5};
};

#endif  // QAG4SIMULATIONDISTORTIONS_H

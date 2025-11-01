// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONDISTORTIONS_H
#define QAG4SIMULATIONDISTORTIONS_H

#include <fun4all/SubsysReco.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeed.h>

#include <math.h>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class TrkrClusterContainer;
class SvtxTrack;
class ActsGeometry;
class PHG4TpcGeomContainer;

class QAG4SimulationDistortions : public SubsysReco
{
 public:
  QAG4SimulationDistortions(const std::string& name = "QAG4SimulationDistortions");

  ~QAG4SimulationDistortions() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  //! track map name
  void set_trackmap_name( const std::string& value )
  { m_trackmapname = value; }

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

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack* track);
  bool checkTrack(SvtxTrack* track);
  bool checkTPOTResidual(SvtxTrack* track);
  SvtxTrackMap* m_trackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TpcGeomContainer *m_tpcGeom = nullptr;

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  int m_event = 0;
  float m_tanAlpha = NAN;
  float m_tanBeta = NAN;
  float m_drphi = NAN;
  float m_dz = NAN;
  float m_clusR = NAN;
  float m_clusPhi = NAN;
  float m_clusZ = NAN;
  float m_clusEta = NAN;
  float m_stateR = NAN;
  float m_statePhi = NAN;
  float m_stateZ = NAN;
  float m_stateEta = NAN;
  float m_stateRPhiErr = NAN;
  float m_stateZErr = NAN;
  float m_clusRPhiErr = NAN;
  float m_clusZErr = NAN;

  int m_layer = -1;
  float m_statePt = NAN;
  float m_statePz = NAN;
  float m_trackPt = NAN;
  float m_trackdEdx = NAN;
  int m_charge = -10;
  int m_crossing = -10;

  TrkrDefs::cluskey m_cluskey = TrkrDefs::CLUSKEYMAX;

  /// disable distortion correction
  bool m_disable_module_edge_corr = false;
  bool m_disable_static_corr = false;
  bool m_disable_average_corr = false;
  bool m_disable_fluctuation_corr = false;

  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_states = 0;
  int m_accepted_states = 0;
  //@}

  bool m_useMicromegas = true;
  double m_minPT = 0.5;

};

#endif  // QAG4SIMULATIONDISTORTIONS_H

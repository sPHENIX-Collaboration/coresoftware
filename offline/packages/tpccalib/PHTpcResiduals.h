#ifndef TRACKRECO_PHTPCRESIDUALS_H
#define TRACKRECO_PHTPCRESIDUALS_H

#include "TpcSpaceChargeMatrixContainer.h"  // for TpcSpaceChargeMa...

#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <Acts/Definitions/Algebra.hpp>  // for Vector3

#include <cmath>
#include <memory>
#include <string>   // for basic_string
#include <utility>  // for pair

class ActsGeometry;
class SvtxTrackState;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrClusterContainer;

/**
 * This class takes preliminary fits from PHActsTrkFitter to the
 * silicon + MM clusters and calculates the residuals in the TPC
 * from that track fit. The TPC state has to be explicitly determined
 * here since the Acts::DirectNavigator does not visit the TPC states
 */
class PHTpcResiduals : public SubsysReco
{
 public:
  PHTpcResiduals(const std::string &name = "PHTpcResiduals");
  ~PHTpcResiduals() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  ///@name Option for setting distortion correction calculation limits
  //@{
  void setMaxTrackAlpha(float maxTAlpha)
  {
    m_maxTAlpha = maxTAlpha;
  }

  void setMaxTrackBeta(float maxTBeta)
  {
    m_maxTBeta = maxTBeta;
  }

  void setMaxTrackResidualDrphi(float maxResidualDrphi)
  {
    m_maxResidualDrphi = maxResidualDrphi;
  }

  void setMaxTrackResidualDz(float maxResidualDz)
  {
    m_maxResidualDz = maxResidualDz;
  }
  //@}

  void setMinRPhiErr(float minRPhiErr)
  {
    m_minRPhiErr = minRPhiErr;
  }

  void setMinZErr(float minZErr)
  {
    m_minZErr = minZErr;
  }

  /// track min pT
  void setMinPt(double value)
  {
    m_minPt = value;
  }

  /// track crossing
  void requireCrossing(bool flag = true)
  {
    m_requireCrossing = flag;
  }

  /// Grid dimensions
  void setGridDimensions(const int phiBins, const int rBins, const int zBins);

  /// set to true to store evaluation histograms and ntuples
  void setSavehistograms(bool) {}

  /// output file name for evaluation histograms
  void setHistogramOutputfile(const std::string &) {}

  /// output file name for storing the space charge reconstruction matrices
  void setOutputfile(const std::string &outputfile)
  {
    m_outputfile = outputfile;
  }

  /// require micromegas to be present when extrapolating tracks to the TPC
  void setUseMicromegas(bool value)
  {
    m_useMicromegas = value;
  }

  void disableModuleEdgeCorr() { m_disable_module_edge_corr = true; }
  void disableStaticCorr() { m_disable_static_corr = true; }
  void disableAverageCorr() { m_disable_average_corr = true; }
  void disableFluctuationCorr() { m_disable_fluctuation_corr = true; }

  /// modify track map name
  void setTrackMapName( const std::string& value )
  { m_trackmapname = value; }

 private:

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  int processTracks(PHCompositeNode *topNode);

  bool checkTrack(SvtxTrack *track) const;
  bool checkTPOTResidual(SvtxTrack *track) const;
  void processTrack(SvtxTrack *track);

  /// Gets distortion cell for identifying bins in TPC
  int getCell(const Acts::Vector3 &loc);

  /// Node information for Acts tracking geometry and silicon+MM
  /// track fit
  std::string m_trackmapname = "SvtxSiliconMMTrackMap";
  SvtxTrackMap *m_trackMap = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  float m_maxTAlpha = 0.6;
  float m_maxResidualDrphi = 0.5;  // cm
  float m_maxTBeta = 1.5;
  float m_maxResidualDz = 0.5;  // cm

  float m_minRPhiErr = 0.005;  // 0.005cm -- 50um
  float m_minZErr = 0.01;      // 0.01cm -- 100um

  static constexpr float m_phiMin = 0;
  static constexpr float m_phiMax = 2. * M_PI;

  static constexpr float m_rMin = 20;  // cm
  static constexpr float m_rMax = 78;  // cm

  static constexpr int m_minClusCount = 10;

  static constexpr float m_layerMin = 7;
  static constexpr float m_layerMax = 55;

  /// Tpc geometry
  static constexpr unsigned int m_nLayersTpc = 48;
  float m_zMin = 0;  // cm
  float m_zMax = 0;  // cm

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

  // TODO: check if needed
  int m_event = 0;

  /// require micromegas to be present when extrapolating tracks to the TPC
  bool m_useMicromegas = true;

  /// minimum pT required for track to be considered in residuals calculation (GeV/c)
  double m_minPt = 0.5;

  /// require track crossing zero
  bool m_requireCrossing = false;

  /// disable distortion correction
  bool m_disable_module_edge_corr = false;
  bool m_disable_static_corr = false;
  bool m_disable_average_corr = false;
  bool m_disable_fluctuation_corr = false;

  /// output file
  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_clusters = 0;
  int m_accepted_clusters = 0;
  //@}
};

#endif

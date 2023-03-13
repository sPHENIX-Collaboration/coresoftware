#ifndef TRACKRECO_PHTPCRESIDUALS_H
#define TRACKRECO_PHTPCRESIDUALS_H

#include <fun4all/SubsysReco.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

#include <Acts/Utilities/Result.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <memory>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TpcSpaceChargeMatrixContainer;
class TrkrCluster;
class TrkrClusterContainer;

class TFile;
class TH1;
class TH2;
class TTree;

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

  /// Option for setting distortion correction calculation limits
  void setMaxTrackAlpha(float maxTAlpha) 
    { m_maxTAlpha = maxTAlpha;}
  void setMaxTrackBeta(float maxTBeta)
    { m_maxTBeta = maxTBeta; }
  void setMaxTrackResidualDrphi(float maxResidualDrphi) 
    { m_maxResidualDrphi = maxResidualDrphi;}
  
  void setMaxTrackResidualDz(float maxResidualDz)
    { m_maxResidualDz = maxResidualDz; }
  
  void setGridDimensions(const int phiBins, const int rBins, const int zBins);

  /// set to true to store evaluation histograms and ntuples
  void setSavehistograms( bool value ) { m_savehistograms = value; }
    
  /// output file name for evaluation histograms
  void setHistogramOutputfile(const std::string &outputfile) {m_histogramfilename = outputfile;}

  /// output file name for storing the space charge reconstruction matrices
  void setOutputfile(const std::string &outputfile) {m_outputfile = outputfile;}

  /// require micromegas to be present when extrapolating tracks to the TPC
  void setUseMicromegas( bool value )
  { m_useMicromegas = value; }

  /// cluster version
  /* Note: this could be retrived automatically using dynamic casts from TrkrCluster objects */
  void setClusterVersion(int value) { m_cluster_version = value; }

  private:

  using BoundTrackParam = 
    const Acts::BoundTrackParameters;
  
  /// pairs path length and track parameters
  using BoundTrackParamPair = std::pair<float,BoundTrackParam>;

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*, short int crossing) const;

  int processTracks(PHCompositeNode *topNode);

  bool checkTrack(SvtxTrack* track) const;
  void processTrack(SvtxTrack* track);

  /// fill track state from bound track parameters
  void addTrackState( SvtxTrack* track, float pathlength, const Acts::BoundTrackParameters& params );
  
  /** \brief 
   * Propagates the silicon+MM track fit to the surface on which
   * an available source link in the TPC exists, added from the stub
   * matching propagation
   * returns the path lenght and the resulting parameters
   */
  BoundTrackParamPair propagateTrackState( const Acts::BoundTrackParameters& params, const Surface& surf ) const;

  /// Gets distortion cell for identifying bins in TPC
  int getCell(const Acts::Vector3& loc);

  /// create histograms
  void makeHistograms();
  
  Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track) const;

  /// actis transformation
  ActsTransformations m_transformer;
  
  /// Node information for Acts tracking geometry and silicon+MM
  /// track fit
  SvtxTrackMap *m_trackMap = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;

  // crossing z correction
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  
  // distortion corrections
  TpcDistortionCorrectionContainer* m_dcc_static = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_average = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_fluctuation = nullptr;

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;
  
  float m_maxTAlpha = 0.6;
  float m_maxResidualDrphi = 0.5; // cm
  float m_maxTBeta = 1.5;
  float m_maxResidualDz = 0.5; // cm

  static constexpr float m_phiMin = 0;
  static constexpr float m_phiMax = 2. * M_PI;

  static constexpr float m_rMin = 20; // cm
  static constexpr float m_rMax = 78; // cm

  static constexpr int m_minClusCount = 10;

  /// Tpc geometry
  static constexpr unsigned int m_nLayersTpc = 48;
  static constexpr float m_zMin = -105.5; // cm
  static constexpr float m_zMax = 105.5;  // cm
 
  /// cluster error parametrisation
  ClusterErrorPara m_cluster_error_parametrization;
  
  /// cluster version
  int m_cluster_version = 4;

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;
  
  // TODO: check if needed
  int m_event = 0;
  
  /// require micromegas to be present when extrapolating tracks to the TPC
  bool m_useMicromegas = true;

  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  /// running track crossing id
  short int m_crossing = 0;
  
  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_clusters = 0;
  int m_accepted_clusters = 0;
  //@}

  /// Output root histograms
  bool m_savehistograms = false;
  TH2 *h_rphiResid = nullptr;
  TH2 *h_zResid = nullptr;
  TH2 *h_etaResidLayer = nullptr;
  TH2 *h_zResidLayer = nullptr;
  TH2 *h_etaResid = nullptr;
  TH1 *h_index = nullptr;
  TH2 *h_alpha = nullptr;
  TH2 *h_beta = nullptr;
  
  //@name additional histograms that copy the per-cell data used to extract the distortions
  //@{
  using TH1_map_t = std::map<int,TH1*>;
  using TH2_map_t = std::map<int,TH2*>;
  
  TH1_map_t h_drphi;
  TH1_map_t h_dz;
  TH2_map_t h_drphi_alpha;
  TH2_map_t h_dz_beta;
  //@}
  
  TTree *residTup = nullptr;

  /// delta rphi vs layer number
  TH2 *h_deltarphi_layer = nullptr;

  /// delta z vs layer number
  TH2 *h_deltaz_layer = nullptr;

  std::string m_histogramfilename = "PHTpcResiduals.root";
  std::unique_ptr<TFile> m_histogramfile = nullptr;

  /// For diagnostics
  double tanAlpha = 0;
  double tanBeta = 0;
  double drphi = 0;
  double dz = 0;
  double clusR = 0;
  double clusPhi = 0;
  double clusZ = 0;
  double statePhi = 0;
  double stateZ = 0;
  double stateRPhiErr = 0;
  double stateZErr = 0;
  double clusRPhiErr = 0;
  double clusZErr = 0;
  double stateR = 0;
  TrkrDefs::cluskey cluskey = 0;
};

#endif


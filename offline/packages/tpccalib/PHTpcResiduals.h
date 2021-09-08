#ifndef TRACKRECO_PHTPCRESIDUALS_H
#define TRACKRECO_PHTPCRESIDUALS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Utilities/Result.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <ActsExamples/EventData/Track.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrClusterContainer;
class TpcSpaceChargeMatrixContainer;
class TrkrCluster;

namespace ActsExamples
{
  class TrkrClusterSourceLink;
}

#include <memory>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;

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
  
  void setGridDimensions(const int phiBins, const int rBins,
			 const int zBins);

  /// set to true to store evaluation histograms and ntuples
  void setSavehistograms( bool value ) { m_savehistograms = value; }
    
  /// output file name for storing the space charge reconstruction matrices
  void setOutputfile(const std::string &outputfile) {m_outputfile = outputfile;}
  
 private:

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  int processTracks(PHCompositeNode *topNode);

  bool checkTrack(SvtxTrack* track);
  void processTrack(SvtxTrack* track);

  /// Calculates TPC residuals given an Acts::Propagation result to
  /// a TPC surface
  void calculateTpcResiduals(const Acts::BoundTrackParameters& params,
			     const TrkrCluster* cluster);

  /// Propagates the silicon+MM track fit to the surface on which
  /// an available source link in the TPC exists, added from the stub
  /// matching propagation
  BoundTrackParamPtrResult propagateTrackState(
  const Acts::BoundTrackParameters& params, 
		     const SourceLink& sl);

  /// Gets distortion cell for identifying bins in TPC
  int getCell(const Acts::Vector3D& loc);
  
  void makeHistograms();
  SourceLink makeSourceLink(TrkrCluster* cluster);
  Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track);
  Surface getSurface(TrkrDefs::cluskey cluskey,
		     TrkrDefs::subsurfkey);
      
  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey);
  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey);
  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey);
  Acts::Vector3D getVertex(SvtxTrack *track);

  /// Node information for Acts tracking geometry and silicon+MM
  /// track fit
  SvtxVertexMap *m_vertexMap = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;

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

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;
  
  // TODO: check if needed
  int m_event = 0;
  
  /// Counter for number of bad propagations from propagateTrackState()
  int m_nBadProps = 0;

  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

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
  TTree *residTup = nullptr;
  std::string m_histogramfilename = "PHTpcResiduals.root";
  std::unique_ptr<TFile> m_histogramfile = nullptr;

  /// For diagnostics
  double tanAlpha = 0, tanBeta = 0, drphi = 0, dz = 0, clusR = 0, clusPhi = 0, 
    clusZ = 0, statePhi = 0, stateZ = 0, stateRPhiErr = 0, stateZErr = 0, 
    clusRPhiErr = 0, clusZErr = 0, stateR = 0;
  // int cell = 0, ir = 0, iz = 0, iphi = 0;
  unsigned int cluskey = 0;
};

#endif


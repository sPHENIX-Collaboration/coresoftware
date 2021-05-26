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
class TrkrClusterContainer;
class SvtxTrack;
class SvtxVertexMap;
class TrkrCluster;
class SvtxTrackMap;

namespace ActsExamples
{
  class TrkrClusterSourceLink;
}

#include <memory>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
  ~PHTpcResiduals() override;

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

  /// Option for outputting some basic cluster-track 
  /// distortion histograms
  void setOutputfile(std::string outputfile) 
    {m_outputfile = outputfile;}
  
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
  int getCell(const int iz, const int ir, const int iphi);

  /// Calculates distortion matrices and updates pointer to node
  void calculateDistortions(PHCompositeNode *topNode);
  
  void makeHistograms();
  TH3* createHistogram(TH3* hin, const TString& name);
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

  const float m_phiMin = 0;
  const float m_phiMax = 2. * M_PI;

  const float m_rMin = 20; // cm
  const float m_rMax = 78; // cm

  const int m_minClusCount = 10;

  /// Tpc geometry
  const unsigned int m_nLayersTpc = 48;
  const float m_zMin = -105.5; // cm
  const float m_zMax = 105.5;  // cm

  /// These are grid sizes given by the distortion model
  int m_zBins = 80;
  int m_phiBins = 36;
  int m_rBins = 16;
  int m_totalBins = m_zBins * m_phiBins * m_rBins;

  /// Number of dimensions for the residuals below
  const int m_nCoord = 3;

  /// Vectors for collecting TPC residual information
  std::vector<Acts::SymMatrix3D> m_lhs;
  std::vector<Acts::Vector3D> m_rhs;
  std::vector<int> m_clusterCount;
  
  int m_event = 0;
  
  /// Counter for number of bad propagations from propagateTrackState()
  int m_nBadProps = 0;

  /// Output root histograms
  std::string m_outputfile = "PHTpcResidualsDistortionCorrections.root";
  TH2 *h_rphiResid = nullptr;
  TH2 *h_zResid = nullptr;
  TH2 *h_etaResidLayer = nullptr;
  TH2 *h_zResidLayer = nullptr;
  TH2 *h_etaResid = nullptr;
  TH1 *h_index = nullptr;
  TH2 *h_alpha = nullptr;
  TH2 *h_beta = nullptr;
  TTree *residTup = nullptr;
  TFile *m_outputFile = nullptr;

  /// For diagnostics
  double tanAlpha = 0, tanBeta = 0, drphi = 0, dz = 0, clusR = 0, clusPhi = 0, 
    clusZ = 0, statePhi = 0, stateZ = 0, stateRPhiErr = 0, stateZErr = 0, 
    clusRPhiErr = 0, clusZErr = 0, stateR = 0;
  int cell = 0, ir = 0, iz = 0, iphi = 0;
  unsigned int cluskey = 0;
};

#endif


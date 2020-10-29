#ifndef TRACKRECO_PHTPCRESIDUALS_H
#define TRACKRECO_PHTPCRESIDUALS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include "ActsTrack.h"
#include "ActsTrackingGeometry.h"

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Utilities/Result.hpp>

class PHCompositeNode;

#include <memory>
#include <map>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>


using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;

using DistortionMap = std::map<const int, 
                               std::pair<unsigned int, 
                                         const double>>;

/**
 * A struct containing the distortion map correction information
 * to be put on the node tree for other modules
 */
struct DistortionCorrections
{
  DistortionMap *m_distortionMap = nullptr;
  DistortionMap *m_distortionMapErr = nullptr;
  
  /// For definitions see equivalent definitions in PHTpcResiduals
  int m_zBins = 50;
  int m_phiBins = 72;
  int m_rBins = 48;
  const int m_totalBins = m_zBins * m_phiBins * m_rBins;
  int m_nCoord = 3;
};



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
  ~PHTpcResiduals();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  /// Option for setting distortion correction calculation limits
  void setMaxTrackAlpha(float maxTAlpha) { m_maxTAlpha = maxTAlpha;}
  void setMaxTrackResidual(float maxResidual) 
    { m_maxResidual = maxResidual;}
  
  /// Option for outputting some basic cluster-track 
  /// distortion histograms
  void setOutputRoot(bool outputRoot) {m_outputRoot = outputRoot;}
  
 private:

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  int processTracks(PHCompositeNode *topNode);
  void processTrack(ActsTrack& track);

  /// Calculates TPC residuals given an Acts::Propagation result to
  /// a TPC surface
  void calculateTpcResiduals(const Acts::BoundTrackParameters& params,
			     const SourceLink& sl);

  /// Propagates the silicon+MM track fit to the surface on which
  /// an available source link in the TPC exists, added from the stub
  /// matching propagation
  BoundTrackParamPtrResult propagateTrackState(
                     const ActsExamples::TrackParameters& params, 
		     const SourceLink& sl);

  /// Gets distortion cell for identifying bins in TPC
  int getCell(const int actsLayer, const Acts::Vector3D& loc);
  int getCell(const int iz, const int ir, const int iphi);

  /// Calculates distortion matrices and updates pointer to node
  void calculateDistortions(PHCompositeNode *topNode);
  
  void makeHistograms();
  
  /// Node information for Acts tracking geometry and silicon+MM
  /// track fit
  std::map<unsigned int, ActsTrack> *m_actsProtoTracks = nullptr;
  ActsTrackingGeometry *m_tGeometry;
 
  DistortionMap m_distortions;
  DistortionMap m_distortionErrs;

  float m_maxTAlpha = 0.6;
  float m_maxResidual = 5 * Acts::UnitConstants::cm;

  /// Tpc geometry
  const unsigned int m_nLayersTpc = 48;
  const float m_zMin = -2120 / 2.; // mm
  const float m_zMax = 2120 / 2.; // mm

  /// These are grid sizes given by the distortion model
  const int m_zBins = 50;
  const int m_phiBins = 72;
  const int m_rBins = 48;
  const int m_totalBins = m_zBins * m_phiBins * m_rBins;

  /// Number of dimensions for the residuals below
  const int m_nCoord = 3;

  /// Vectors for collecting TPC residual information
  std::vector<Acts::SymMatrix3D> m_lhs;
  std::vector<Acts::Vector3D> m_rhs;
  std::vector<int> m_clusterCount;
  
  int m_event = 0;
  
  /// Counter for number of bad propagations from propagateTrackState()
  int m_nBadProps = 0;

  /// Container for distortion corrections for node tree
  DistortionCorrections *m_distortionCorrections;

  /// Output root histograms
  int m_outputRoot = false;
  TFile *outfile;
  TH2 *h_rphiResid;
  TH2 *h_zResid;
  TH2 *h_etaResidLayer;
  TH2 *h_zResidLayer;
  TH2 *h_etaResid;

};

#endif


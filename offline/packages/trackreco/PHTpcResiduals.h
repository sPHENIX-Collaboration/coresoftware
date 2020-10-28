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
  int ResetEvent(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void setMaxTrackAlpha(float maxTAlpha) { m_maxTAlpha = maxTAlpha;}
  void setMaxTrackResidual(float maxResidual) 
    { m_maxResidual = maxResidual;}
  void setOutputRoot(bool outputRoot) {m_outputRoot = outputRoot;}
  
 private:

  int getNodes(PHCompositeNode *topNode);

  int processTracks(PHCompositeNode *topNode);
  void processTrack(ActsTrack& track);

  void calculateTpcResiduals(const Acts::BoundTrackParameters& params,
			     const SourceLink& sl);

  BoundTrackParamPtrResult propagateTrackState(
                     const ActsExamples::TrackParameters& params, 
		     const SourceLink& sl);

  int getCell(const int actsLayer, const Acts::Vector3D& loc);
  int getCell(const int iz, const int ir, const int iphi);
  void calculateDistortions(PHCompositeNode *topNode);

  std::map<unsigned int, ActsTrack> *m_actsProtoTracks = nullptr;
  ActsTrackingGeometry *m_tGeometry;
 
  float m_maxTAlpha = 0.6;
  float m_maxResidual = 5;

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

  int m_nBadProps = 0;

  int m_outputRoot = false;
  
  TFile *outfile;
  TH2 *h_rphiResid;
  TH2 *h_zResid;
  TH2 *h_etaResidLayer;
  TH2 *h_zResidLayer;
  TH2 *h_etaResid;

};

#endif


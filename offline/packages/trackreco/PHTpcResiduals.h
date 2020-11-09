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
#include <TH1.h>
#include <TH2.h>


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
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

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

  bool checkTrack(ActsTrack& track);
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
  std::map<unsigned int, ActsTrack> *m_actsProtoTracks;
  ActsTrackingGeometry *m_tGeometry = nullptr;
 
  float m_maxTAlpha = 0.6;
  float m_maxResidualDrphi = 5.; // mm
  float m_maxTBeta = 1.5;
  float m_maxResidualDz = 5.; // mm

  const float m_phiMin = 0;
  const float m_phiMax = 2. * M_PI;

  const float m_rMin = 200; // mm
  const float m_rMax = 780; // mm

  const int m_minClusCount = 10;

  /// Tpc geometry
  const unsigned int m_nLayersTpc = 48;
  const float m_zMin = -1050.; // mm
  const float m_zMax = 1050.;  // mm

  /// These are grid sizes given by the distortion model
  int m_zBins = 50;
  int m_phiBins = 72;
  int m_rBins = 48;
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

};

#endif


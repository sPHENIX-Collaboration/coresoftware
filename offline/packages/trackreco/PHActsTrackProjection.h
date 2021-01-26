#ifndef TRACKRECO_PHACTSTRACKPROJECTION_H
#define TRACKRECO_PHACTSTRACKPROJECTION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>

#include "ActsTrack.h"
#include "ActsTrackingGeometry.h"

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class SvtxTrackMap;
class SvtxTrack;

#include <memory>
#include <map>
#include <string>

using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

using FitParameters = Acts::SingleBoundTrackParameters<Acts::SinglyCharged>;


/**
 * This class takes final fitted tracks from the Acts track fitting
 * and projects them out to cylinders with radius at the same radius
 * as the three calorimeters. Cluster matching is performed with the
 * projections and the SvtxTrack object is updated.
 */

class PHActsTrackProjection : public SubsysReco
{

 public:
  PHActsTrackProjection(const std::string& name 
			= "PHActsTrackProjection");
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:
  
  int getNodes(PHCompositeNode *topNode);
  int projectTracks(PHCompositeNode *topNode, int caloLayer);

  /// Propagate the fitted track parameters to a surface with Acts
  BoundTrackParamPtrResult propagateTrack(
	const FitParameters& params, 
	const SurfacePtr &targetSurf);

  /// Set the particular calo nodes depending on which layer
  int setCaloContainerNodes(PHCompositeNode *topNode,
			     const int caloLayer);
  
  /// Make Acts::CylinderSurface objects corresponding to the calos
  int makeCaloSurfacePtrs(PHCompositeNode *topNode);

  /// Update the SvtxTrack object with the track-cluster match
  void updateSvtxTrack(const Acts::BoundTrackParameters& params,
		       const unsigned int trackKey,
		       const int caloLayer);

  /// Get 3x3 and 5x5 tower sums matched to a track
  void getSquareTowerEnergies(int phiBin, int etaBin,
			      double& energy3x3,
			      double& energy5x5);

  /// Get the cluster values for a particular matched track
  void getClusterProperties(double phi, double eta,
			    double& minIndex, double& minDphi,
			    double& minDeta, double& minE);

  double deltaPhi(const double& phi);

  /// Objects containing the Acts track fit results
  ActsTrackingGeometry *m_tGeometry = nullptr;
  std::map<const unsigned int, Trajectory> *m_actsFitResults;
  SvtxTrackMap *m_trackMap = nullptr;
  
  /// Objects to hold calorimeter information. There are 
  /// only 3 calo layers
  const static int m_nCaloLayers = 3;
  std::vector<std::string> m_caloNames;
  std::vector<SvtxTrack::CAL_LAYER> m_caloTypes;
  std::map<std::string, SurfacePtr> m_caloSurfaces;
  
  RawTowerGeomContainer *m_towerGeomContainer = nullptr;
  RawTowerContainer *m_towerContainer = nullptr;
  RawClusterContainer *m_clusterContainer = nullptr;

  bool m_useCemcPosRecalib = false;

  int m_event = 0;
};

#endif

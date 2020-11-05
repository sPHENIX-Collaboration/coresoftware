#ifndef TRACKRECO_PHACTSTRACKPROJECTION_H
#define TRACKRECO_PHACTSTRACKPROJECTION_H


#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

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

#include <memory>
#include <map>
#include <string>

using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

using FitParameters = Acts::SingleBoundTrackParameters<Acts::SinglyCharged>;

class PHActsTrackProjection : public SubsysReco
{

 public:
  PHActsTrackProjection(const std::string& name);
  
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:
  
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  int projectTracks(PHCompositeNode *topNode, int calLayer);

  /// Propagate the fitted track parameters to a surface
  BoundTrackParamPtrResult propagateTrack(
	const FitParameters& params, 
	const SurfacePtr &targetSurf);

  int setCaloContainerNodes(PHCompositeNode *topNode,
			     const int calLayer);
  int makeCaloSurfacePtrs(PHCompositeNode *topNode);

  void updateSvtxTrack(const Acts::BoundTrackParameters& params,
		       const unsigned int trackKey);

  ActsTrackingGeometry *m_tGeometry = nullptr;
  std::map<const unsigned int, Trajectory> *m_actsFitResults;
  int m_event = 0;

  const static int m_nCaloLayers = 3;
  std::string m_caloNames[m_nCaloLayers];
  std::map<std::string, SurfacePtr> m_caloSurfaces;
  
  RawTowerGeomContainer *m_towerGeomContainer = nullptr;
  RawTowerContainer *m_towerContainer = nullptr;
  RawClusterContainer *m_clusterContainer = nullptr;

};

#endif

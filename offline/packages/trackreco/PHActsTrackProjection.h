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
  int projectTracks(PHCompositeNode *topNode, int caloLayer);

  /// Propagate the fitted track parameters to a surface
  BoundTrackParamPtrResult propagateTrack(
	const FitParameters& params, 
	const SurfacePtr &targetSurf);

  int setCaloContainerNodes(PHCompositeNode *topNode,
			     const int caloLayer);
  int makeCaloSurfacePtrs(PHCompositeNode *topNode);

  void updateSvtxTrack(const Acts::BoundTrackParameters& params,
		       const unsigned int trackKey,
		       const int caloLayer);

  void getSquareTowerEnergies(int phiBin, int etaBin,
			      double& energy3x3,
			      double& energy5x5);

  void getClusterProperties(double phi, double eta,
			    double& minIndex, double& minDphi,
			    double& minDeta, double& minE);

  double deltaPhi(const double& phi);

  ActsTrackingGeometry *m_tGeometry = nullptr;
  std::map<const unsigned int, Trajectory> *m_actsFitResults;
  int m_event = 0;
  SvtxTrackMap *m_trackMap;
  
  const static int m_nCaloLayers = 3;
  std::vector<std::string> m_caloNames;
  std::vector<SvtxTrack::CAL_LAYER> m_caloTypes;
  std::map<std::string, SurfacePtr> m_caloSurfaces;
  
  RawTowerGeomContainer *m_towerGeomContainer = nullptr;
  RawTowerContainer *m_towerContainer = nullptr;
  RawClusterContainer *m_clusterContainer = nullptr;

};

#endif

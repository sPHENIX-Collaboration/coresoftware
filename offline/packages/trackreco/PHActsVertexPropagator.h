#ifndef TRACKRECO_PHACTSVERTEXPROPAGATOR_H
#define TRACKRECO_PHACTSVERTEXPROPAGATOR_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsTrackingGeometry.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxTrack;

using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class PHActsVertexPropagator : public SubsysReco
{

 public: 
  PHActsVertexPropagator(const std::string& name = "PHActsVertexPropagator");
  
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
 
 private:

  int getNodes(PHCompositeNode *topNode);
  void setTrackVertexTo0();
  BoundTrackParamPtrResult propagateTrack(const Acts::BoundTrackParameters& params,
					  const unsigned int vtxid);
  Acts::Vector3 getVertex(const unsigned int vtxid);
  void updateSvtxTrack(SvtxTrack* track, const Acts::BoundTrackParameters& params);
  void setVtxChi2();
  
  ActsTrackingGeometry *m_tGeometry = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  std::map<const unsigned int, Trajectory> *m_trajectories = nullptr;

};

#endif 

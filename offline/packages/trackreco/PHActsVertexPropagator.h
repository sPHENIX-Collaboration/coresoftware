#ifndef TRACKRECO_PHACTSVERTEXPROPAGATOR_H
#define TRACKRECO_PHACTSVERTEXPROPAGATOR_H

#include "ActsPropagator.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxTrack;

class PHActsVertexPropagator : public SubsysReco
{
 public:
  using BoundTrackParam =
      const Acts::BoundTrackParameters;
  using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
  using SurfacePtr = std::shared_ptr<const Acts::Surface>;
  using Trajectory = ActsExamples::Trajectories;

  PHActsVertexPropagator(const std::string &name = "PHActsVertexPropagator");

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void fieldMap(std::string &fieldmap) { m_fieldMap = fieldmap; }

 private:
  int getNodes(PHCompositeNode *topNode);
  ActsPropagator::BTPPairResult
  propagateTrack(const Acts::BoundTrackParameters &params,
                 const unsigned int vtxid);
  Acts::Vector3 getVertex(const unsigned int vtxid);
  void updateSvtxTrack(SvtxTrack *track,
                       const Acts::BoundTrackParameters &params);
  void setVtxChi2();

  ActsGeometry *m_tGeometry = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  std::map<const unsigned int, Trajectory> *m_trajectories = nullptr;
  std::string m_fieldMap = "";
};

#endif

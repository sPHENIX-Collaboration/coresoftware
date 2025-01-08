
#include "PHActsTrackPropagator.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsAborter.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

//____________________________________________________________________________..
PHActsTrackPropagator::PHActsTrackPropagator(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHActsTrackPropagator::~PHActsTrackPropagator()
{
}

//____________________________________________________________________________..
int PHActsTrackPropagator::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsTrackPropagator::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);

  return ret;
}

//____________________________________________________________________________..
int PHActsTrackPropagator::process_event(PHCompositeNode *)
{
  ActsPropagator prop(m_tGeometry);
  for (auto &[key, track] : *m_trackMap)
  {
    auto params = prop.makeTrackParams(track, m_vertexMap);
    if(!params.ok())
      {
	continue;
      }
    auto result = propagateTrack(params.value());
    if (result.ok())
    {
      addTrackState(result, track);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrackPropagator::addTrackState(
    BoundTrackParamResult &result,
    SvtxTrack *svtxTrack)
{
  float pathlength = result.value().first / Acts::UnitConstants::cm;
  auto params = result.value().second;

  SvtxTrackState_v1 out(pathlength);

  auto projectionPos = params.position(m_tGeometry->geometry().getGeoContext());
  const auto momentum = params.momentum();
  out.set_x(projectionPos.x() / Acts::UnitConstants::cm);
  out.set_y(projectionPos.y() / Acts::UnitConstants::cm);
  out.set_z(projectionPos.z() / Acts::UnitConstants::cm);
  out.set_px(momentum.x());
  out.set_py(momentum.y());
  out.set_pz(momentum.z());

  if (Verbosity() > 1)
  {
    std::cout << "Adding track state for layer " << m_sphenixLayer
              << " with path length " << pathlength << " with position "
              << projectionPos.transpose() << std::endl;
  }

  ActsTransformations transformer;
  const auto globalCov = transformer.rotateActsCovToSvtxTrack(params);
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      out.set_error(i, j, globalCov(i, j));
    }
  }

  svtxTrack->insert_state(&out);
}

PHActsTrackPropagator::BoundTrackParamResult
PHActsTrackPropagator::propagateTrack(const Acts::BoundTrackParameters &params)
{
  ActsPropagator propagator(m_tGeometry);
  propagator.verbosity(Verbosity());

  return propagator.propagateTrack(params, m_sphenixLayer);
}

//____________________________________________________________________________..
int PHActsTrackPropagator::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PHActsTrackPropagator::Print(const std::string &what) const
{
  std::cout << "PHActsTrackPropagator:: " << what << std::endl;
}

int PHActsTrackPropagator::getNodes(PHCompositeNode *topNode)
{
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No vertex map on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(
      topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsTrackingGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "No SvtxTrackMap on node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

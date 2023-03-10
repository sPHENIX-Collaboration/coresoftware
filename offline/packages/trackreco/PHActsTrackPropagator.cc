
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
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
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
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  ret = checkLayer();

  return ret;
}

//____________________________________________________________________________..
int PHActsTrackPropagator::process_event(PHCompositeNode *)
{
  for (auto &[key, track] : *m_trackMap)
  {
    const auto params = makeTrackParams(track);

    auto result = propagateTrack(params);
    if (!std::isnan(result.first))
    {
      addTrackState(result, track);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrackPropagator::addTrackState(
    const BoundTrackParamResult &result,
    SvtxTrack *svtxTrack)
{
  float pathlength = result.first;
  auto params = result.second;

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
  if (Verbosity() > 1)
  {
    std::cout << "Propagating final track fit with momentum: "
              << params.momentum() << " and position "
              << params.position(m_tGeometry->geometry().getGeoContext())
              << std::endl
              << "track fit phi/eta "
              << atan2(params.momentum()(1),
                       params.momentum()(0))
              << " and "
              << atanh(params.momentum()(2) / params.momentum().norm())
              << std::endl;
  }

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;

  auto field = m_tGeometry->geometry().magField;
  auto trackingGeometry = m_tGeometry->geometry().tGeometry;
  Stepper stepper(field);
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);

  Propagator propagator(stepper, navigator);

  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  if (Verbosity() > 3)
  {
    logLevel = Acts::Logging::VERBOSE;
  }

  auto logger = Acts::getDefaultLogger("PHActsTrackPropagator",
                                       logLevel);
  using Actors = Acts::ActionList<>;
  using Aborters = Acts::AbortList<ActsAborter>;

  Acts::PropagatorOptions<Actors, Aborters> options(
      m_tGeometry->geometry().getGeoContext(),
      m_tGeometry->geometry().magFieldContext,
      Acts::LoggerWrapper{*logger});

  options.abortList.get<ActsAborter>().abortlayer = m_actslayer;
  options.abortList.get<ActsAborter>().abortvolume = m_actsvolume;

  auto result = propagator.propagate(params, options);

  if (result.ok())
  {
    auto params = *result.value().endParameters;
    double pathlength = result.value().pathLength / Acts::UnitConstants::cm;
    return std::make_pair(pathlength, params);
  }

  return std::make_pair(NAN, Acts::BoundTrackParameters(nullptr, Acts::BoundVector::Zero(), -1));
}

Acts::BoundTrackParameters
PHActsTrackPropagator::makeTrackParams(SvtxTrack *track)
{
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());

  auto actsVertex = getVertex(track);
  auto perigee =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);
  auto actsFourPos =
      Acts::Vector4(track->get_x() * Acts::UnitConstants::cm,
                    track->get_y() * Acts::UnitConstants::cm,
                    track->get_z() * Acts::UnitConstants::cm,
                    10 * Acts::UnitConstants::ns);

  ActsTransformations transformer;

  Acts::BoundSymMatrix cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsTrackFittingAlgorithm::TrackParameters::create(perigee,
                                                            m_tGeometry->geometry().getGeoContext(),
                                                            actsFourPos, momentum,
                                                            track->get_charge() / track->get_p(),
                                                            cov)
      .value();
}
Acts::Vector3 PHActsTrackPropagator::getVertex(SvtxTrack *track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex *svtxVertex = m_vertexMap->get(vertexId);
  Acts::Vector3 vertex = Acts::Vector3::Zero();
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
    vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
    vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
  }

  return vertex;
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

int PHActsTrackPropagator::checkLayer()
{
  /*
   * Acts geometry is defined in terms of volumes and layers. Within a volume
   * layers always begin at 2 and iterate in 2s, i.e. the MVTX is defined as a
   * volume and the 3 layers are identifiable as 2, 4, and 6.
   * So we convert the sPHENIX layer number here to the Acts volume and
   * layer number that can interpret where to navigate to in the propagation.
   * The only exception is the TPOT, which is interpreted as a single layer.
   */

  /// mvtx
  if (m_sphenixLayer < 3)
  {
    m_actsvolume = 10;
    m_actslayer = (m_sphenixLayer + 1) * 2;
  }

  /// intt
  else if (m_sphenixLayer < 7)
  {
    m_actsvolume = 12;
    m_actslayer = ((m_sphenixLayer - 3) + 1) * 2;
  }

  /// tpc
  else if (m_sphenixLayer < 55)
  {
    m_actsvolume = 14;
    m_actslayer = ((m_sphenixLayer - 7) + 1) * 2;
  }
  /// tpot only has one layer in Acts geometry
  else
  {
    m_actsvolume = 16;
    m_actslayer = 2;
  }

  /// Test to make sure the found volume and layer exist
  auto tgeometry = m_tGeometry->geometry().tGeometry;

  bool foundlayer = false;
  bool foundvolume = false;

  tgeometry->visitSurfaces([&](const Acts::Surface *srf)
                           {
    if (srf != nullptr) {
      auto layer = srf->geometryId().layer();
      auto volume = srf->geometryId().volume();
      if(volume == m_actsvolume and layer == m_actslayer)
	{
	  foundlayer = true;
	  foundvolume = true;
	  return;
	}
    } });

  if (!foundlayer or !foundvolume)
  {
    std::cout << PHWHERE
              << "Could not identify Acts volume and layer to propagate to. Can't continue. "
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

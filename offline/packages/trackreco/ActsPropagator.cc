#include "ActsPropagator.h"

#include <trackbase/ActsAborter.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Propagator/detail/VoidPropagatorComponents.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

ActsPropagator::SurfacePtr
ActsPropagator::makeVertexSurface(const SvtxVertex* vertex)
{
  return Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3(vertex->get_x() * Acts::UnitConstants::cm,
                    vertex->get_y() * Acts::UnitConstants::cm,
                    vertex->get_z() * Acts::UnitConstants::cm));
}
ActsPropagator::SurfacePtr
ActsPropagator::makeVertexSurface(const Acts::Vector3& vertex)
{
  return Acts::Surface::makeShared<Acts::PerigeeSurface>(
      vertex * Acts::UnitConstants::cm);
}
ActsPropagator::BoundTrackParamResult
ActsPropagator::makeTrackParams(SvtxTrackState* state,
                                int trackCharge,
                                ActsPropagator::SurfacePtr surf)
{
  Acts::Vector3 momentum(state->get_px(),
                         state->get_py(),
                         state->get_pz());
  Acts::Vector4 actsFourPos = Acts::Vector4(state->get_x() * Acts::UnitConstants::cm,
                                            state->get_y() * Acts::UnitConstants::cm,
                                            state->get_z() * Acts::UnitConstants::cm,
                                            10 * Acts::UnitConstants::ns);

  ActsTransformations transformer;
  Acts::BoundSquareMatrix cov = transformer.rotateSvtxTrackCovToActs(state);

  return ActsTrackFittingAlgorithm::TrackParameters::create(
             surf,
             m_geometry->geometry().getGeoContext(),
             actsFourPos, momentum,
             trackCharge / momentum.norm(),
             cov,
             Acts::ParticleHypothesis::pion());
}
ActsPropagator::BoundTrackParamResult
ActsPropagator::makeTrackParams(SvtxTrack* track,
                                SvtxVertexMap* vertexMap)
{
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());

  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = vertexMap->get(vertexId);
  Acts::Vector3 vertex = Acts::Vector3::Zero();
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
    vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
    vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
  }

  auto perigee =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(vertex);
  auto actsFourPos =
      Acts::Vector4(track->get_x() * Acts::UnitConstants::cm,
                    track->get_y() * Acts::UnitConstants::cm,
                    track->get_z() * Acts::UnitConstants::cm,
                    10 * Acts::UnitConstants::ns);

  ActsTransformations transformer;

  Acts::BoundSquareMatrix cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsTrackFittingAlgorithm::TrackParameters::create(perigee,
                                                            m_geometry->geometry().getGeoContext(),
                                                            actsFourPos, momentum,
                                                            track->get_charge() / track->get_p(),
                                                            cov,
                                                            Acts::ParticleHypothesis::pion());
}

ActsPropagator::BTPPairResult
ActsPropagator::propagateTrack(const Acts::BoundTrackParameters& params,
                               const unsigned int sphenixLayer)
{
  unsigned int actsvolume, actslayer;
  if (!checkLayer(sphenixLayer, actsvolume, actslayer) or !m_geometry)
  {
    return Acts::Result<BoundTrackParamPair>::failure(std::error_code(0, std::generic_category()));
  }

  if (m_verbosity > 1)
  {
    printTrackParams(params);
  }

  auto propagator = makePropagator();

  using Actors = Acts::ActionList<>;
  using Aborters = Acts::AbortList<ActsAborter>;

  Acts::PropagatorOptions<Actors, Aborters> options(
      m_geometry->geometry().getGeoContext(),
      m_geometry->geometry().magFieldContext);

  options.abortList.get<ActsAborter>().abortlayer = actslayer;
  options.abortList.get<ActsAborter>().abortvolume = actsvolume;

  auto result = propagator.propagate(params, options);

  if (result.ok())
  {
    auto finalparams = *result.value().endParameters;
    auto pathlength = result.value().pathLength;
    auto pair = std::make_pair(pathlength, finalparams);

    return Acts::Result<BoundTrackParamPair>::success(pair);
  }

  return result.error();
}

ActsPropagator::BTPPairResult
ActsPropagator::propagateTrack(const Acts::BoundTrackParameters& params,
                               const SurfacePtr& surface)
{
  if (m_verbosity > 1)
  {
    printTrackParams(params);
  }

  auto propagator = makePropagator();

  Acts::PropagatorOptions<> options(m_geometry->geometry().getGeoContext(),
                                    m_geometry->geometry().magFieldContext);

  auto result = propagator.propagate(params, *surface,
                                     options);

  if (result.ok())
  {
    auto finalparams = *result.value().endParameters;
    auto pathlength = result.value().pathLength;
    auto pair = std::make_pair(pathlength, finalparams);

    return Acts::Result<BoundTrackParamPair>::success(pair);
  }

  return result.error();
}

ActsPropagator::BTPPairResult
ActsPropagator::propagateTrackFast(const Acts::BoundTrackParameters& params,
                                   const SurfacePtr& surface)
{
  if (m_verbosity > 1)
  {
    printTrackParams(params);
  }

  auto propagator = makeFastPropagator();

  Acts::PropagatorOptions<> options(m_geometry->geometry().getGeoContext(),
                                    m_geometry->geometry().magFieldContext);

  auto result = propagator.propagate(params, *surface,
                                     options);

  if (result.ok())
  {
    auto finalparams = *result.value().endParameters;
    auto pathlength = result.value().pathLength;
    auto pair = std::make_pair(pathlength, finalparams);

    return Acts::Result<BoundTrackParamPair>::success(pair);
  }

  return result.error();
}

ActsPropagator::FastPropagator ActsPropagator::makeFastPropagator()
{
  auto field = m_geometry->geometry().magField;

  if (m_constField)
  {
    if (m_verbosity > 2)
    {
      std::cout << "Using const field of val " << m_fieldval << std::endl;
    }
    Acts::Vector3 fieldVec(0, 0, m_fieldval);
    field = std::make_shared<Acts::ConstantBField>(fieldVec);
  }

  ActsPropagator::Stepper stepper(field);

  Acts::Logging::Level logLevel = Acts::Logging::FATAL;
  if (m_verbosity > 3)
  {
    logLevel = Acts::Logging::VERBOSE;
  }

  std::shared_ptr<const Acts::Logger> logger = Acts::getDefaultLogger("ActsPropagator", logLevel);

  return ActsPropagator::FastPropagator(stepper, Acts::detail::VoidNavigator(),
                                        logger);
}
ActsPropagator::SphenixPropagator ActsPropagator::makePropagator()
{
  auto field = m_geometry->geometry().magField;

  if (m_constField)
  {
    Acts::Vector3 fieldVec(0, 0, m_fieldval);
    field = std::make_shared<Acts::ConstantBField>(fieldVec);
  }

  auto trackingGeometry = m_geometry->geometry().tGeometry;
  Stepper stepper(field, m_overstepLimit);
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);

  Acts::Logging::Level logLevel = Acts::Logging::FATAL;
  if (m_verbosity > 3)
  {
    logLevel = Acts::Logging::VERBOSE;
  }

  std::shared_ptr<const Acts::Logger> logger = Acts::getDefaultLogger("ActsPropagator", logLevel);
  return SphenixPropagator(stepper, navigator, logger);
}

bool ActsPropagator::checkLayer(const unsigned int& sphenixlayer,
                                unsigned int& actsvolume,
                                unsigned int& actslayer)
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
  if (sphenixlayer < 3)
  {
    actsvolume = 10;
    actslayer = (sphenixlayer + 1) * 2;
  }

  /// intt
  else if (sphenixlayer < 7)
  {
    actsvolume = 12;
    actslayer = ((sphenixlayer - 3) + 1) * 2;
  }

  /// tpc
  else if (sphenixlayer < 55)
  {
    actsvolume = 14;
    actslayer = ((sphenixlayer - 7) + 1) * 2;
  }
  /// tpot only has one layer in Acts geometry
  else
  {
    actsvolume = 16;
    actslayer = 2;
  }

  /// Test to make sure the found volume and layer exist
  auto tgeometry = m_geometry->geometry().tGeometry;

  bool foundlayer = false;
  bool foundvolume = false;

  tgeometry->visitSurfaces([&](const Acts::Surface* srf)
                           {
    if (srf != nullptr) {
      auto layer = srf->geometryId().layer();
      auto volume = srf->geometryId().volume();
      if(volume == actsvolume and layer == actslayer)
	{
	  foundlayer = true;
	  foundvolume = true;
	  return;
	}
    } });

  if (!foundlayer or !foundvolume)
  {
    std::cout << "trackreco::ActsPropagator::checkLayer Could not identify Acts volume and layer to propagate to. Can't continue. "
              << std::endl;

    return false;
  }

  return true;
}

void ActsPropagator::printTrackParams(const Acts::BoundTrackParameters& params)
{
  std::cout << "Propagating final track fit with momentum: "
            << params.momentum() << " and position "
            << params.position(m_geometry->geometry().getGeoContext())
            << std::endl
            << "track fit phi/eta "
            << atan2(params.momentum()(1),
                     params.momentum()(0))
            << " and "
            << atanh(params.momentum()(2) / params.momentum().norm())
            << std::endl;
}


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


#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/ActsAborter.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Propagator/Navigator.hpp>

//____________________________________________________________________________..
PHActsTrackPropagator::PHActsTrackPropagator(const std::string &name):
 SubsysReco(name)
{

}

//____________________________________________________________________________..
PHActsTrackPropagator::~PHActsTrackPropagator()
{

}

//____________________________________________________________________________..
int PHActsTrackPropagator::Init(PHCompositeNode*)
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
int PHActsTrackPropagator::process_event(PHCompositeNode*)
{
 
  for(auto &[key, track] : *m_trackMap)
    {
      const auto params = makeTrackParams(track); 

      auto result = propagateTrack(params);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


PHActsTrackPropagator::BoundTrackParamResult
PHActsTrackPropagator::propagateTrack(const Acts::BoundTrackParameters& params)
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
  using Aborter = ActsAborter;
  using Aborters = Acts::AbortList<ActsAborter>;
  using abortlist = Acts::AbortList<>;
  Acts::PropagatorOptions<Actors, Aborters> options(
     m_tGeometry->geometry().getGeoContext(),
     m_tGeometry->geometry().magFieldContext,
     Acts::LoggerWrapper{*logger});

  auto result = propagator.propagate(params, options);
  if(result.ok())
    {
      return Acts::Result<BoundTrackParam>::success(std::move((*result).endParameters.value()));
    }

  return result.error();

}

Acts::BoundTrackParameters
PHActsTrackPropagator::makeTrackParams(SvtxTrack* track)
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
     cov).value();
}
Acts::Vector3 PHActsTrackPropagator::getVertex(SvtxTrack* track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
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
int PHActsTrackPropagator::End(PHCompositeNode*)
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

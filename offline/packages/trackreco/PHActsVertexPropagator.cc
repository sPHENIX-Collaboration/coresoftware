#include "PHActsVertexPropagator.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

PHActsVertexPropagator::PHActsVertexPropagator(const std::string& name)
  : SubsysReco(name)
{}


int PHActsVertexPropagator::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexPropagator::InitRun(PHCompositeNode* topNode)
{
  int returnval = getNodes(topNode);
  return returnval;
}
int PHActsVertexPropagator::process_event(PHCompositeNode*)
{
  std::vector<unsigned int> deletedKeys;
  for(const auto& [trackKey, trajectory] : *m_trajectories)
    {
      auto svtxTrack = m_trackMap->get(trackKey);
      if(!svtxTrack)
	{
	  /// Key was removed by the track cleaner, remove it from
	  /// the trajectory list too
	  deletedKeys.push_back(trackKey);
	  continue;
	}
      if(Verbosity() > 2)
	{ svtxTrack->identify(); }
      const auto &[trackTips, mj] = trajectory.trajectory();

      if(trackTips.size() > 1 and Verbosity() > 0)
	{ 
	  std::cout << PHWHERE 
		    << "More than 1 track tip per track. Should never happen..."
		    << std::endl;
	}
      
      for(const auto& trackTip : trackTips)
	{
	  const auto& boundParams = trajectory.trackParameters(trackTip);
	 
	  auto propresult = propagateTrack(boundParams, svtxTrack->get_vertex_id());
	  if(propresult.ok())
	    {
	  
	      auto paramsAtVertex = std::move(**propresult);
	      updateSvtxTrack(svtxTrack,paramsAtVertex);
	    }
	}
    }

  /// Erase the trajectories that were removed from the track cleaner
  for(auto& key : deletedKeys)
    {
      m_trajectories->erase(key);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsVertexPropagator::updateSvtxTrack(SvtxTrack* track, 
					     const Acts::BoundTrackParameters& params)
{
  auto position = params.position(m_tGeometry->geoContext);
  
  if(Verbosity() > 2)
    {
      std::cout << "Updating position track parameters from " << track->get_x()
		<< ", " << track->get_y() << ", " << track->get_z() << " to " 
		<< position.transpose() / 10.
		<< std::endl;
    }

  track->set_x(position(0) / Acts::UnitConstants::cm);
  track->set_y(position(1) / Acts::UnitConstants::cm);
  track->set_z(position(2) / Acts::UnitConstants::cm);

}

BoundTrackParamPtrResult PHActsVertexPropagator::propagateTrack(
		         const Acts::BoundTrackParameters& params,
			 const unsigned int vtxid)
{
  
  /// create perigee surface
  auto actsVertex = getVertex(vtxid);
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);

  return std::visit(
      [params, perigee, this]
      (auto && inputField) ->BoundTrackParamPtrResult {
	using InputMagneticField = 
	  typename std::decay_t<decltype(inputField)>::element_type;
	using MagneticField      = Acts::SharedBField<InputMagneticField>;
	using Stepper            = Acts::EigenStepper<MagneticField>;
	using Propagator         = Acts::Propagator<Stepper>;
	
	MagneticField field(inputField);
	Stepper stepper(field);
	Propagator propagator(stepper);
	
	Acts::Logging::Level logLevel = Acts::Logging::FATAL;
	if(Verbosity() > 3)
	  { logLevel = Acts::Logging::VERBOSE; }
	
	auto logger = Acts::getDefaultLogger("PHActsVertexPropagator", 
					     logLevel);
	
	Acts::PropagatorOptions<> options(m_tGeometry->geoContext,
					  m_tGeometry->magFieldContext,
					  Acts::LoggerWrapper{*logger});
	
	auto result = propagator.propagate(params, *perigee, 
					   options);
	if(result.ok())
	  { return std::move((*result).endParameters); }
	
	return result.error();
      },
      m_tGeometry->magField);
}

Acts::Vector3D PHActsVertexPropagator::getVertex(const unsigned int vtxid)
{
  auto svtxVertex = m_vertexMap->get(vtxid);
  return Acts::Vector3D(svtxVertex->get_x() * Acts::UnitConstants::cm,
			svtxVertex->get_y() * Acts::UnitConstants::cm,
			svtxVertex->get_z() * Acts::UnitConstants::cm);
}

int PHActsVertexPropagator::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsVertexPropagator::getNodes(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "Acts tracking geometry not on node tree, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");

  if(!m_trajectories)
    {
      std::cout << PHWHERE << "No acts trajectories on node tree, exiting. " 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "No svtx vertex map, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "No svtx track map, exiting. " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

#include "PHActsVertexPropagator.h"
#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

PHActsVertexPropagator::PHActsVertexPropagator(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
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

  if(m_vertexMap->size() == 0)
    { setTrackVertexTo0(); }

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
      
      const auto& trackTips = trajectory.tips();
     
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

  ActsTransformations rotater;
  rotater.setVerbosity(Verbosity());
  if(params.covariance())
    {
      auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params, m_tGeometry->geoContext);
      
      /// Update covariance
      for(int i = 0; i < 3; i++) {
	for(int j = 0; j < 3; j++) {
	  track->set_error(i,j, rotatedCov(i,j));
	  if(i < 2 and j < 2)
	    { track->set_acts_covariance(i,j, params.covariance().value()(i,j)); }
	}
      }
    }

  updateTrackDCA(track);

}

void PHActsVertexPropagator::updateTrackDCA(SvtxTrack* track)
{
  Acts::Vector3 pos(track->get_x(),
		    track->get_y(),
		    track->get_z());
  Acts::Vector3 mom(track->get_px(),
		    track->get_py(),
		    track->get_pz());

  auto vtxid = track->get_vertex_id();
  auto svtxVertex = m_vertexMap->get(vtxid);
  Acts::Vector3 vertex(svtxVertex->get_x(),
		       svtxVertex->get_y(),
		       svtxVertex->get_z());

  pos -= vertex;

  Acts::ActsSymMatrix<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i, j) = track->get_error(i, j);
	} 
    }
  
  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi = atan2(r(1), r(0));
  
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  
  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * pos;
  Acts::ActsSymMatrix<3> rotCov = rot * posCov * rot_T;

  const auto dca3Dxy = pos_R(0);
  const auto dca3Dz = pos_R(2);
  const auto dca3DxyCov = rotCov(0,0);
  const auto dca3DzCov = rotCov(2,2);

  track->set_dca3d_xy(dca3Dxy);
  track->set_dca3d_z(dca3Dz);
  track->set_dca3d_xy_error(sqrt(dca3DxyCov));
  track->set_dca3d_z_error(sqrt(dca3DzCov));

}
BoundTrackParamPtrResult PHActsVertexPropagator::propagateTrack(
		         const Acts::BoundTrackParameters& params,
			 const unsigned int vtxid)
{
  
  /// create perigee surface
  auto actsVertex = getVertex(vtxid);
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  
  Stepper stepper(m_tGeometry->magField);
  Acts::Navigator::Config cfg{m_tGeometry->tGeometry};
  Acts::Navigator navigator(cfg);
  Propagator propagator(stepper, navigator);
  
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

}

Acts::Vector3 PHActsVertexPropagator::getVertex(const unsigned int vtxid)
{
  auto svtxVertex = m_vertexMap->get(vtxid);
  return Acts::Vector3(svtxVertex->get_x() * Acts::UnitConstants::cm,
		       svtxVertex->get_y() * Acts::UnitConstants::cm,
		       svtxVertex->get_z() * Acts::UnitConstants::cm);
}

void PHActsVertexPropagator::setTrackVertexTo0()
{
  /// If we found no vertices in the event, propagate the tracks to 0,0,0
  auto vertex = std::make_unique<SvtxVertex_v1>();
  vertex->set_chisq(0.);
  vertex->set_ndof(0);
  vertex->set_t0(0);
  vertex->set_id(0);
  vertex->set_x(0);
  vertex->set_y(0);
  vertex->set_z(0);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      vertex->set_error(i,j, 20.);
    }
  }

  m_vertexMap->insert(vertex.release());

  for(auto& [key, track] : *m_trackMap)
    {
      track->set_vertex_id(0);
    }

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

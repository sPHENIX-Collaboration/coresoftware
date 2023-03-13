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

	  auto result = propagateTrack(boundParams, svtxTrack->get_vertex_id());
	  if(result.ok())
	    {
	      updateSvtxTrack(svtxTrack, result.value());
	    }
	  else
	    {
	      svtxTrack->identify();
	    }
	}
    }
  
  setVtxChi2();

  /// Erase the trajectories that were removed from the track cleaner
  for(auto& key : deletedKeys)
    {
      m_trajectories->erase(key);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsVertexPropagator::setVtxChi2()
{

  for(const auto& [vtxid, vtx] : *m_vertexMap)
    {
      float xvtx = vtx->get_x();
      float yvtx = vtx->get_y();
      float zvtx = vtx->get_z();
      float xchisqsum = 0;
      float ychisqsum = 0;
      float zchisqsum = 0;
      
      for(auto trackiter = vtx->begin_tracks(); trackiter != vtx->end_tracks();
	  ++trackiter)
	{
	  SvtxTrack* track = m_trackMap->get(*trackiter);
	  if(!track) { continue; }

	  float trkx = track->get_x();
	  float trky = track->get_y();
	  float trkz = track->get_z();
	  float trkcovx = track->get_error(0,0);
	  float trkcovy = track->get_error(1,1);
	  float trkcovz = track->get_error(2,2);
	  
	  xchisqsum += pow(trkx - xvtx, 2) / trkcovx;
	  ychisqsum += pow(trky - yvtx, 2) / trkcovy;
	  zchisqsum += pow(trkz - zvtx, 2) / trkcovz;

	}

      /// independent chisq sum additively
      vtx->set_chisq(xchisqsum + ychisqsum + zchisqsum);
      /// Each track contributes independently to x,y,z, so the total
      /// ndf is total tracks * 3 minus 1*3 for each independent x,y,z
      vtx->set_ndof(vtx->size_tracks() * 3 - 3);
    }
 
}

void PHActsVertexPropagator::updateSvtxTrack(SvtxTrack* track, 
					     const Acts::BoundTrackParameters& params)
{
  auto position = params.position(m_tGeometry->geometry().getGeoContext());
  
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
      auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params);
    
      /// Update covariance
      for(int i = 0; i < 3; i++) {
	for(int j = 0; j < 3; j++) {
	  track->set_error(i,j, rotatedCov(i,j));
	}
      }
    }
}

BoundTrackParamResult PHActsVertexPropagator::propagateTrack(
		           const Acts::BoundTrackParameters& params,
			   const unsigned int vtxid)
{
  
  /// create perigee surface
  auto actsVertex = getVertex(vtxid);
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  
  Stepper stepper(m_tGeometry->geometry().magField);
  Acts::Navigator::Config cfg{m_tGeometry->geometry().tGeometry};
  Acts::Navigator navigator(cfg);
  Propagator propagator(stepper, navigator);
  
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  if(Verbosity() > 3)
    { logLevel = Acts::Logging::VERBOSE; }
  
  auto logger = Acts::getDefaultLogger("PHActsVertexPropagator", 
				       logLevel);

  Acts::PropagatorOptions<> options(m_tGeometry->geometry().getGeoContext(),
				    m_tGeometry->geometry().magFieldContext,
				    Acts::LoggerWrapper{*logger});
  
  auto result = propagator.propagate(params, *perigee, 
				     options);
  
  if(result.ok())
    { 
      return Acts::Result<BoundTrackParam>::success(std::move((*result).endParameters.value()));
      return params;
    }

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
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
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

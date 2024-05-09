#include "PHActsToSvtxTracks.h"
#include <trackbase_historic/ActsTransformations.h>

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHTimer.h>

#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/Utilities/Units.hpp>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include <TMatrixDSym.h>


PHActsToSvtxTracks::PHActsToSvtxTracks(const std::string &name)
  : SubsysReco(name)
  , m_actsFitResults(nullptr)
{
  Verbosity(0);
}

int PHActsToSvtxTracks::End(PHCompositeNode *topNode)
{

  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsToSvtxTracks" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsToSvtxTracks::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsToSvtxTracks::InitRun(PHCompositeNode *topNode)
{
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsToSvtxTracks::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Start process_event in PHActsToSvtxTracks" << std::endl;
  }

  for(auto& [trackKey, trajectory] : *m_actsFitResults)
    {
      createSvtxTrack(trackKey, trajectory);
    }
  
  if(Verbosity() > 0)
    std::cout << "Finished PHActsToSvtxTracks process_event"
	      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsToSvtxTracks::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsToSvtxTracks::createSvtxTrack(const unsigned int trackKey,
					 Trajectory traj)
{
  /// Try to find the track first to update it. Otherwise, 
  /// if it is a new map, create the track
  SvtxTrack *svtxTrack = m_svtxTrackMap->find(trackKey)->second;

  const auto vertexId = svtxTrack->get_vertex_id();

  auto svtxVertex = m_svtxVertexMap->find(vertexId)->second;

  Acts::Vector3D vertex(svtxVertex->get_x() * Acts::UnitConstants::cm,
		      svtxVertex->get_y() * Acts::UnitConstants::cm,
		      svtxVertex->get_z() * Acts::UnitConstants::cm);
 
  const auto &[trackTips, mj] = traj.trajectory();
  
  /// For a track from the Acts KF, it has only one track tip
  const auto& trackTip = trackTips.front();
  const auto& params = traj.trackParameters(trackTip);
  
  svtxTrack->clear_states();
  svtxTrack->clear_cluster_keys();
  
  float pathlength = 0.0;
  SvtxTrackState_v1 out( pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  svtxTrack->insert_state(&out);   
  
  auto trajState = 
    Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
  
  svtxTrack->set_x(params.position(m_tGeometry->getGeoContext())(0)
	       / Acts::UnitConstants::cm);
  svtxTrack->set_y(params.position(m_tGeometry->getGeoContext())(1)
	       / Acts::UnitConstants::cm);
  svtxTrack->set_z(params.position(m_tGeometry->getGeoContext())(2)
	       / Acts::UnitConstants::cm);
  
  svtxTrack->set_px(params.momentum()(0));
  svtxTrack->set_py(params.momentum()(1));
  svtxTrack->set_pz(params.momentum()(2));
  
  svtxTrack->set_charge(params.charge());
  svtxTrack->set_chisq(trajState.chi2Sum);
  svtxTrack->set_ndf(trajState.NDF);

  float dca3Dxy = NAN;
  float dca3Dz = NAN;
  float dca3DxyCov = NAN;
  float dca3DzCov = NAN;
  
  auto rotater = std::make_unique<ActsTransformations>();
  rotater->setVerbosity(Verbosity());
  
  if(params.covariance())
    {
      auto rotatedCov = 
	rotater->rotateActsCovToSvtxTrack(params,
					  m_tGeometry->getGeoContext());
      
      for(int i = 0; i < 6; i++)
	for(int j = 0; j < 6; j++)
	  svtxTrack->set_error(i, j, rotatedCov(i,j));
      
      rotater->calculateDCA(params, vertex, rotatedCov,
			    m_tGeometry->getGeoContext(),
			    dca3Dxy, dca3Dz,
			    dca3DxyCov, dca3DzCov);
      
    }
       
  svtxTrack->set_dca3d_xy(dca3Dxy / Acts::UnitConstants::cm);
  svtxTrack->set_dca3d_z(dca3Dz / Acts::UnitConstants::cm);

  /// Units already considered in rotation
  svtxTrack->set_dca3d_xy_error(sqrt(dca3DxyCov));
  svtxTrack->set_dca3d_z_error(sqrt(dca3DzCov));
  
  rotater->fillSvtxTrackStates(traj, trackTip, svtxTrack,
			       m_tGeometry->getGeoContext());
  return;
}

int PHActsToSvtxTracks::createNodes(PHCompositeNode *topNode)
{
  
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
  }
  
  PHCompositeNode *svtxNode = 
    dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, 
						    m_svtxMapName.c_str());
  if(!m_svtxTrackMap)
    {
      m_svtxTrackMap = new SvtxTrackMap_v1;
      PHIODataNode<PHObject> *trackNode = 
	new PHIODataNode<PHObject>(m_svtxTrackMap,m_svtxMapName.c_str(),
				   "PHObject");
      dstNode->addNode(trackNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;;
}

int PHActsToSvtxTracks::getNodes(PHCompositeNode *topNode)
{
 

  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!m_svtxTrackMap)
  {
    std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_svtxVertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_svtxVertexMap)
    {
      std::cout << PHWHERE << "SvtxVertexMap not found on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");

  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");

  if(!m_actsFitResults)
    {
      std::cout << PHWHERE << "ActsFitResults not on node tree, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


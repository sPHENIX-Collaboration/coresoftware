
#include "PHTrackClusterAssociator.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>

#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <phgeom/PHGeomUtility.h>


//____________________________________________________________________________..
PHTrackClusterAssociator::PHTrackClusterAssociator(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
PHTrackClusterAssociator::~PHTrackClusterAssociator()
{
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) 
    { return ret; }

  ret = createNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    { return ret; }

  return ret;
}


//____________________________________________________________________________..
int PHTrackClusterAssociator::process_event(PHCompositeNode *topNode)
{

  if(!m_calosAvailable) { return Fun4AllReturnCodes::EVENT_OK; }

  for(int layer = 0; layer < m_nCaloLayers; layer++)
    {
      if(Verbosity() > 1)
	{ 
	  std::cout << "Processing calo layer " 
		  << m_caloNames.at(layer) << std::endl;
	}

      int ret = matchTracks(topNode, layer);
      if(ret != Fun4AllReturnCodes::EVENT_OK)
	{ return Fun4AllReturnCodes::ABORTEVENT; }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackClusterAssociator::matchTracks(PHCompositeNode* topNode, 
					  const int caloLayer)
{

  if (getCaloNodes(topNode, caloLayer) 
      != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }

  for(const auto& [key, track] : *m_trackMap)
    {
      float caloRadius = m_caloRadii.at(m_caloTypes.at(caloLayer));
      const SvtxTrackState* state = track->get_state(caloRadius);
      const float statex = state->get_x();
      const float statey = state->get_y();
      const float statez = state->get_z();
      const float statephi = atan2(statey, statex);
      const float stateeta = asinh(statez / 
			     sqrt(statex * statex 
				  + statey * statey));


     


    }





}
int PHTrackClusterAssociator::getNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!m_trackMap)
    {
      std::cout << "No track map, can't continue" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHTrackClusterAssociator::getCaloNodes(PHCompositeNode *topNode,
					   const int caloLayer)
{
  std::string towerGeoNodeName = "TOWERGEOM_" + m_caloNames.at(caloLayer);
  std::string towerNodeName    = "TOWER_CALIB_" + m_caloNames.at(caloLayer);
  std::string clusterNodeName  = "CLUSTER_" + m_caloNames.at(caloLayer);

  m_towerGeomContainer = findNode::getClass<RawTowerGeomContainer>
    (topNode, towerGeoNodeName.c_str());

  m_towerContainer = findNode::getClass<RawTowerContainer>
    (topNode, towerNodeName.c_str());
  
  m_clusterContainer = findNode::getClass<RawClusterContainer>
    (topNode, clusterNodeName.c_str());

  if(m_useCemcPosRecalib and 
     m_caloNames.at(caloLayer).compare("CEMC") == 0) 
    {
      std::string nodeName = "CLUSTER_POS_COR_" + m_caloNames.at(caloLayer);
      m_clusterContainer = findNode::getClass<RawClusterContainer>
	(topNode, nodeName.c_str());
    }
  
  if(!m_towerGeomContainer or !m_towerContainer or !m_clusterContainer)
    {
      std::cout << PHWHERE 
		<< "Calo geometry and/or cluster container not found on node tree. Track-Calo cluster map won't be filled."
		<< std::endl;
      m_calosAvailable = false;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackClusterAssociator::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsSiliconSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_trackClusterMap = findNode::getClass<TrackClusterMap>(topNode, "TrackCaloClusterMap");
  if(!m_trackClusterMap)
    {
      m_trackClusterMap = new TrackClusterMap;
       PHDataNode<TrackClusterMap> *newnode 
	= new PHDataNode<TrackClusterMap>(m_trackClusterMap, "TrackCaloClusterMap");
      svtxNode->addNode(newnode);
    }

  return Fun4AllReturnCodes::EVENT_OK;

} 


//____________________________________________________________________________..
int PHTrackClusterAssociator::ResetEvent(PHCompositeNode*)
{

  m_trackClusterMap->clear();
 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


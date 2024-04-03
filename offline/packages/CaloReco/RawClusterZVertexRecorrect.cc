#include "RawClusterZVertexRecorrect.h"

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>                         // for pair

RawClusterZVertexRecorrect::RawClusterZVertexRecorrect(const std::string &name)
  : SubsysReco(std::string("RawClusterZVertexRecorrect_") + name)
  , _det_name("CEMC")  // not tested for hcal yet
  , m_UseTowerInfo(0)
  , m_UseBbcZVtx(false)
{

}

int RawClusterZVertexRecorrect::InitRun(PHCompositeNode *topNode)
{
  if (!topNode && Verbosity())
  {
    std::cout << "RawClusZVtxRecorrect::InitRun :  NO TOPNODE" << std::endl;
  }

  //  CreateNodeTree(topNode);
  m_calrecoUtilInstance.LoadProfile();
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterZVertexRecorrect::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "RawClusZVtxRecorrect::Processing a NEW EVENT" << std::endl;
  }

  std::string rawClusNodeName = "CLUSTER_" + _det_name;
  if (m_UseTowerInfo)
  {
    rawClusNodeName = "CLUSTERINFO_" + _det_name;
  }

  RawClusterContainer *rawclusters = findNode::getClass<RawClusterContainer>(topNode, rawClusNodeName.c_str());
  if (!rawclusters)
  {
    std::cout << "No " << _det_name << " Cluster Container found while in RawClusterZVertexRecorrect, can't proceed!!!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get Vertex
  //float vx = 0;
  //  float vy = 0;
  float vz = 0;


  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    
  if (vertexmap && m_UseTowerInfo < 0)
    //  if (vertexmap)
    {
      if (!vertexmap->empty())
	{
	  GlobalVertex *vertex = (vertexmap->begin()->second);
	  //	  vx = vertex->get_x();
	  // vy = vertex->get_y();
	  vz = vertex->get_z();
	}
    }
  
  
  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
 
  if (m_UseBbcZVtx && mbdmap && m_UseTowerInfo < 2)
    {
      //      std::cout << " in mbdmap ccpi0 " << std::endl;
      
      MbdVertex *bvertex = nullptr;
      for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
           mbditer != mbdmap->end();
           ++mbditer)
	{

	  bvertex = mbditer->second;
	}
      //      MbdVertex *bvertex = (mbdmap->begin()->second);
      if (!bvertex) 
      {
	return Fun4AllReturnCodes::ABORTEVENT;
      }
      vz = bvertex->get_z();
    }

  /////////////////////////////

  // loop over the clusters
  RawClusterContainer::ConstRange begin_end = rawclusters->getClusters();
  RawClusterContainer::ConstIterator iter;

  for (iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    //    RawClusterDefs::keytype key = iter->first;
    RawCluster *cluster = iter->second;

    float clus_energy = cluster->get_energy();
    float clus_savz = cluster->get_z();
    float clus_chi2 = cluster->get_chi2();

    m_calrecoUtilInstance.ShowerDepthCorrZVertex(cluster, vz);
    m_calrecoUtilInstance.ProbCorrsZVertex(cluster, vz);


    /*
    // for if it is desired in future to make new copied node instead of 
    // changing clusterNode "in Situ"
  
    RawCluster *recalibcluster = dynamic_cast<RawCluster *>(cluster->CloneMe());
    assert(recalibcluster);
    //    if (m_UseTowerInfo)
    //  std::cout << "and here" << std::endl;
    recalibcluster->set_energy(clus_energy / eclus_recalib_val);
    recalibcluster->set_ecore(cluster->get_ecore() / ecore_recalib_val);
    _recalib_clusters->AddCluster(recalibcluster);
    */

    if (Verbosity() && clus_energy > 18.0)
    {
      std::cout << "Input,out eclus cluster energies: " << clus_energy 
		<< " " << cluster->get_energy() << std::endl;
      std::cout << "Input, out  cluster z:" << clus_savz 
		<< " " << cluster->get_z() << std::endl;

      std::cout << "Input, out  cluster ch2:" << clus_chi2 
		<< " " << cluster->get_chi2() << std::endl;
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;
}


/*

  //keeping this for if we want to make a new node like CLUSTER_POS_CORR
void RawClusterZVertexRecorrect::CreateNodeTree(PHCompositeNode *topNode)
{

  // Check that it is there
  if (!dstNode)
  {
    std::cout << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in RawClusterZVertexRecorrect::CreateNodeTree");
  }

  // Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _det_name));

  // Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode(_det_name);
    dstNode->addNode(cemcNode);
  }


  // Check to see if the cluster recalib node is on the nodetree
  _recalib_clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_RECALIB_" + _det_name);
  std::string ClusterCorrNodeName = "CLUSTER_POS_COR_" + _det_name;
  ;

  // If not, make it and add it to the _det_name subnode
  if (!_recalib_clusters)
  {
    _recalib_clusters = new RawClusterContainer();
    if (m_UseTowerInfo)
    {
      ClusterCorrNodeName = "CLUSTERINFO_POS_COR_" + _det_name;
    }

    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_recalib_clusters, ClusterCorrNodeName.c_str(), "PHObject");
    cemcNode->addNode(clusterNode);
  }


}

*/

int RawClusterZVertexRecorrect::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

/*!
 * \file DSTClusterPruning.cc
 * \author Alex Patton <aopatton@mit.edu>
 */

#include "DSTClusterPruning.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterv5.h>
//include new cluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
// #include <trackbase/TrkrClusterContainerv5.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed_v1.h>
//include SvtxTrackMap_v1 to write data to it
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v4.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>
#include <utility> //for pair

#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TNtuple.h>
//_____________________________________________________________________

//_____________________________________________________________________
DSTClusterPruning::DSTClusterPruning( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int DSTClusterPruning::Init(PHCompositeNode* topNode )
{
  std::cout << "Writer Init start" << std::endl;
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTClusterPruning::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "DSTClusterPruning::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if( !svtxNode )
  {
    // create
    std::cout << "DSTTrackArrayReader::Init - SVTX node missing - creating" << std::endl;
    svtxNode = new PHCompositeNode( "SVTX" );
    dstNode->addNode(svtxNode);
  }

  auto trkrNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TRKR"));
  if( !trkrNode )
  {
    // create
    std::cout << "DSTTrackArrayReader::Init - TRKR node missing - creating" << std::endl;
    trkrNode = new PHCompositeNode( "TRKR" );
    dstNode->addNode(trkrNode);
  }

  //make new cluster container
  auto clsNode = findNode::getClass<TrkrClusterContainer>(trkrNode, "ReducedClusterContainer");
  if (!clsNode)
  {
    auto newClusterNode = new PHIODataNode<PHObject>(new TrkrClusterContainerv4, "ReducedClusterContainer", "PHObject");
    trkrNode->addNode(newClusterNode);
  }

                              
  std::cout << "Writer Init end" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTClusterPruning::InitRun(PHCompositeNode* topNode)
{ 
  return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTClusterPruning::process_event(PHCompositeNode* topNode)
{
  //make topNode run in Init
  //Init(topNode);
  // load nodes
  std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << std::endl;
  std::cout << "DSTClusterPruning::process_event" << std::endl;
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) {
    return res;
  }
  std::cout << "Return codes  start" << Fun4AllReturnCodes::EVENT_OK << std::endl;


  if( m_reduced_cluster_map) {
    m_reduced_cluster_map->Reset();
  }
 
  prune_clusters();
  
  std::cout << "Return codes end" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTClusterPruning::End(PHCompositeNode* topNode)
{
  std::cout << "In the end" << std::endl;
  // tcl->Print();
  // tcl->Write();
  // fcl->Write();
  // fcl->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTClusterPruning::load_nodes( PHCompositeNode* topNode )
{


  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  //look for reduced cluster
  m_reduced_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "ReducedClusterContainer");
  
  //m_reduced_track_map = findNode::getClass<SvtxTrackMap>(topNode, "ReducedTrackContainer");

  return Fun4AllReturnCodes::EVENT_OK;

}



//_____________________________________________________________________
void DSTClusterPruning::prune_clusters(){
//use this to create object that looks through both tracks and clusters and saves into new object
//make sure clusters exist
//std::cout << "start of check" << "\n";
//if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
//make sure tracks exist
if( !( m_track_map && m_cluster_map && m_reduced_cluster_map) ) {
    return;
}

  for( const auto& trackpair:*m_track_map )
  {
    std::cout << "start of loop" << "\n";


    //unsigned int key = trackpair.first;
    const auto track = trackpair.second;

    TrackSeed* TPCSeed = track->get_tpc_seed();
    TrackSeed* SiliconSeed = track->get_silicon_seed();

  
    if(!TPCSeed){
      std::cout << "TPCSeed does not exist \n";
      
    }else{
    std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
    TPCSeed->identify();
    for( auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if(!m_reduced_cluster_map->findCluster(cluster_key)){
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_reduced_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        

      }

    }

  }

  
  if(!SiliconSeed){
    std::cout << "SiliconSeed does not exist \n";
    }else{
    std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
  SiliconSeed->identify();
    for( auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if(!m_reduced_cluster_map->findCluster(cluster_key)){
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_reduced_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        
        //m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, cluster);
      }
    
      
    }
  }
   
    
      
      std::cout << "end of loop" << "\n";
    }




}

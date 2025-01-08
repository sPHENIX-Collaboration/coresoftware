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
  auto clsNode = findNode::getClass<TrkrClusterContainer>(trkrNode, "TRKR_CLUSTER_SEED");
  if (!clsNode)
  {
    auto newClusterNode = new PHIODataNode<PHObject>(new TrkrClusterContainerv4, "TRKR_CLUSTER_SEED", "PHObject");
    trkrNode->addNode(newClusterNode);
  }


                              
  std::cout << "Writer Init end" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


//_____________________________________________________________________
//int DSTClusterPruning::InitRun(PHCompositeNode*  /*topNode*/)

//{ 
//  return Fun4AllReturnCodes::EVENT_OK; }
  

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
  //print_clusters();
  //if( m_cluster_map) {
  //  m_cluster_map->Reset();
  //}
  //fill_clusters();
  /*
  if( m_reduced_cluster_map) {
    m_reduced_cluster_map->Reset();
  }
  */
  std::cout << "Return codes end" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
//int DSTClusterPruning::End(PHCompositeNode*  /*topNode*/)
//{
  //std::cout << "In the end" << std::endl;
  // tcl->Print();
  // tcl->Write();
  // fcl->Write();
  // fcl->Close();

  //return Fun4AllReturnCodes::EVENT_OK;
//}

//_____________________________________________________________________
int DSTClusterPruning::load_nodes( PHCompositeNode* topNode )
{


  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  //look for reduced cluster
  m_reduced_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  
  //m_reduced_track_map = findNode::getClass<SvtxTrackMap>(topNode, "ReducedTrackContainer");
  m_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  m_tpc_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  m_silicon_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  //load in the Track Seeds
  //TpcTrackSeedContainer
  //SiliconTrackSeedContainer
  //SvtxTrackSeedContainer

  return Fun4AllReturnCodes::EVENT_OK;

}



//_____________________________________________________________________
void DSTClusterPruning::prune_clusters(){
//use this to create object that looks through both tracks and clusters and saves into new object
//make sure clusters exist
//std::cout << "start of check" << "\n";
//if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
//make sure tracks exist
if( !(m_cluster_map && m_reduced_cluster_map && m_track_seed_container && m_silicon_track_seed_container && m_tpc_track_seed_container) ) {
    return;
}
//std::cout <<"after check " << std::endl;

  //for( auto seed_iter = m_track_seed_container->begin(); seed_iter != m_track_seed_container->end(); ++ seed_iter){
  //const auto& trackseed = *seed_iter;
  for(const auto& trackseed:*m_track_seed_container){
    //std::cout << "Inside trackseed container loop" << std::endl; 
    if(!trackseed){
      continue;
    }
    //std::cout << "after check if trackseed exists" << std::endl;
    //trackseed->identify();
    //trackseed is an SvtxTrackSeed_v2
    unsigned int tpcIndex = trackseed->get_tpc_seed_index();
    unsigned int siliconIndex = trackseed->get_silicon_seed_index();

    TrackSeed* TPCSeed = nullptr;
    if(tpcIndex <= m_tpc_track_seed_container->end()-m_tpc_track_seed_container->begin()){
      TPCSeed = *(m_tpc_track_seed_container->begin()+tpcIndex);
    }
    TrackSeed* SiliconSeed = nullptr;
    if(siliconIndex <= m_silicon_track_seed_container->end()-m_silicon_track_seed_container->begin()){
    SiliconSeed = *(m_silicon_track_seed_container->begin()+siliconIndex);
    }


    if(!TPCSeed){
      std::cout << "TPCSeed does not exist \n";
      
    }else{
    //std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
    //TPCSeed->identify();
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
        //std::cout<< "Cluster RPhi Error: " << cluster->getRPhiError() << std::endl;
        //std::cout<< "Cluster Z Error: " << cluster->getZError() << std::endl;
        //PrintCluster(cluster);
        //cluster->identify();
        // std::cout << "Cluster PhiSize: " << cluster->getPhiSize() << std::endl;
        // std::cout << "Cluster ZSize: " << cluster->getZSize() << std::endl;
        // std::cout << "Cluster RPhiError: " << cluster->getRPhiError() << std::endl;
        // std::cout << "Cluster ZError: " << cluster->getZError() << std::endl;
        // std::cout << "Cluster adc: " << cluster->getAdc() << std::endl;
        // std::cout << "Cluster max adc: " << cluster->getMaxAdc() << std::endl;
        // std::cout << "Cluster subsurfkey: " << cluster->getSubSurfKey() << std::endl;
        // std::cout << "Cluster overlap: " << cluster->getOverlap() << std::endl;
        // std::cout << "Cluster edge: " << cluster->getEdge() << std::endl;
        //cluster->identify();
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_cluster->setPhiError(cluster->getRPhiError());
        m_reduced_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        //PrintCluster(m_cluster);

         // m_cluster->identify();
        // std::cout << "m_Cluster PhiSize: " << m_cluster->getPhiSize() << std::endl;
        // std::cout << "m_Cluster ZSize: " << m_cluster->getZSize() << std::endl;
        // std::cout << "m_Cluster RPhiError: " << m_cluster->getRPhiError() << std::endl;
        // std::cout << "m_Cluster ZError: " << m_cluster->getZError() << std::endl;
        // std::cout << "m_Cluster adc: " << m_cluster->getAdc() << std::endl;
        // std::cout << "m_Cluster max adc: " << m_cluster->getMaxAdc() << std::endl;
        // std::cout << "m_Cluster subsurfkey: " << m_cluster->getSubSurfKey() << std::endl;
        // std::cout << "m_Cluster overlap: " << m_cluster->getOverlap() << std::endl;
        // std::cout << "m_Cluster edge: " << m_cluster->getEdge() << std::endl;
        //std::cout<< "m_Cluster RPhi Error: " << m_cluster->getRPhiError() << std::endl;
        //std::cout<< "m_Cluster Z Error: " << m_cluster->getZError() << std::endl;
        //m_cluster->identify();

      }

    }

  }

  
  if(!SiliconSeed){
    std::cout << "SiliconSeed does not exist \n";
    }else{
    //std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
    //SiliconSeed->identify();
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
        //std::cout<< "Cluster RPhi Error: " << cluster->getRPhiError() << std::endl;
        //std::cout<< "Cluster Z Error: " << cluster->getZError() << std::endl;
        // cluster->identify();
        //PrintCluster(cluster);
        // std::cout << "Cluster PhiSize: " << cluster->getPhiSize() << std::endl;
        // std::cout << "Cluster ZSize: " << cluster->getZSize() << std::endl;
        // std::cout << "Cluster RPhiError: " << cluster->getRPhiError() << std::endl;
        // std::cout << "Cluster ZError: " << cluster->getZError() << std::endl;
        // std::cout << "Cluster adc: " << cluster->getAdc() << std::endl;
        // std::cout << "Cluster max adc: " << cluster->getMaxAdc() << std::endl;
        // std::cout << "Cluster subsurfkey: " << cluster->getSubSurfKey() << std::endl;
        // std::cout << "Cluster overlap: " << cluster->getOverlap() << std::endl;
        // std::cout << "Cluster edge: " << cluster->getEdge() << std::endl;
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_cluster->setPhiError(cluster->getRPhiError());
        m_reduced_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        //PrintCluster(m_cluster);

        // m_cluster->identify();
        // std::cout << "m_Cluster PhiSize: " << m_cluster->getPhiSize() << std::endl;
        // std::cout << "m_Cluster ZSize: " << m_cluster->getZSize() << std::endl;
        // std::cout << "m_Cluster RPhiError: " << m_cluster->getRPhiError() << std::endl;
        // std::cout << "m_Cluster ZError: " << m_cluster->getZError() << std::endl;
        // std::cout << "m_Cluster adc: " << m_cluster->getAdc() << std::endl;
        // std::cout << "m_Cluster max adc: " << m_cluster->getMaxAdc() << std::endl;
        // std::cout << "m_Cluster subsurfkey: " << m_cluster->getSubSurfKey() << std::endl;
        // std::cout << "m_Cluster overlap: " << m_cluster->getOverlap() << std::endl;
        // std::cout << "m_Cluster edge: " << m_cluster->getEdge() << std::endl;

        //std::cout<< "m_Cluster RPhi Error: " << m_cluster->getRPhiError() << std::endl;
        //std::cout<< "m_Cluster Z Error: " << m_cluster->getZError() << std::endl;
        //m_cluster->identify();
        //m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, cluster);
      }
    
      
    }
  }

    //loop over trackseeds and fill clusters
    
    /*for( auto key_iter = trackseed->begin_cluster_keys(); key_iter != trackseed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      std::cout << "cluster key is: " << cluster_key << std::endl;
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

    }*/

  }
  /*
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
    */




}

//fill original clusters
//_____________________________________________________________________
void DSTClusterPruning::fill_clusters(){
//use this to create object that looks through both tracks and clusters and saves into new object
//make sure clusters exist
//std::cout << "start of check" << "\n";
//if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
//make sure tracks exist
if( !(m_cluster_map && m_reduced_cluster_map && m_track_seed_container && m_silicon_track_seed_container && m_tpc_track_seed_container) ) {
    return;
}

  //for( auto seed_iter = m_track_seed_container->begin(); seed_iter != m_track_seed_container->end(); ++ seed_iter){
  //const auto& trackseed = *seed_iter;
  for(const auto& trackseed:*m_track_seed_container){
    //std::cout << "Inside trackseed container loop" << std::endl; 
    if(!trackseed){
      continue;
    }
    //std::cout << "after check if trackseed exists" << std::endl;
    //trackseed->identify();
    //trackseed is an SvtxTrackSeed_v2
    unsigned int tpcIndex = trackseed->get_tpc_seed_index();
    unsigned int siliconIndex = trackseed->get_silicon_seed_index();

    TrackSeed* TPCSeed = nullptr;
    if(tpcIndex <= m_tpc_track_seed_container->end()-m_tpc_track_seed_container->begin()){
      TPCSeed = *(m_tpc_track_seed_container->begin()+tpcIndex);
    }
    TrackSeed* SiliconSeed = nullptr;
    if(siliconIndex <= m_silicon_track_seed_container->end()-m_silicon_track_seed_container->begin()){
    SiliconSeed = *(m_silicon_track_seed_container->begin()+siliconIndex);
    }


    if(!TPCSeed){
      std::cout << "TPCSeed does not exist \n";
      
    }else{
    //std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
    //TPCSeed->identify();
    for( auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_reduced_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if(!m_cluster_map->findCluster(cluster_key)){
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_cluster->setPhiError(cluster->getRPhiError());
        m_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        

      }

    }

  }

  
  if(!SiliconSeed){
    std::cout << "SiliconSeed does not exist \n";
    }else{
    //std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
    //SiliconSeed->identify();
    for( auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_reduced_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if(!m_cluster_map->findCluster(cluster_key)){
        m_cluster = new TrkrClusterv5();
        m_cluster->CopyFrom(cluster);
        m_cluster->setPhiError(cluster->getRPhiError());
        m_cluster_map->addClusterSpecifyKey(cluster_key,m_cluster);
        
        //m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, cluster);
      }
    
      
    }
  }

  }



}




//print clusters to test
void DSTClusterPruning::print_clusters(){

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
      auto reducedcluster = m_reduced_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if( !reducedcluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find reducedcluster for key " << cluster_key << std::endl;
        continue;
      }
      std::cout << "ClusterKey: " << cluster_key << std::endl;
      std::cout << "Cluster map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
      std::cout << "Reduced map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
    

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

      auto reducedcluster = m_reduced_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      if( !reducedcluster )
      {
        std::cout << "DSTClusterPruning::evaluate_tracks - unable to find reducedcluster for key " << cluster_key << std::endl;
        continue;
      }
      std::cout << "ClusterKey: " << cluster_key << std::endl;
      std::cout << "Cluster map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
      std::cout << "Reduced map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
    
      
    }
  }
   
    
      
      std::cout << "end of loop" << "\n";
    }

}

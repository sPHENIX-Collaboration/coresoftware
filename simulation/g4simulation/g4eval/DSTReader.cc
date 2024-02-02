/*!
 * \file DSTReader.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "DSTReader.h"
#include "DSTContainerv3.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <trackbase/InttDefs.h>
#include <intt/InttClusterizer.h>
#include <micromegas/MicromegasDefs.h>
#include <trackbase/MvtxDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
//#include <trackbase/TrkrCluster.h>
//#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterContainerv5.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
//#include <trackbase_historic/TrackSeed_v1.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>
#define __STDC_LIMIT_MACROS

//_____________________________________________________________________
namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    inline const typename T::first_type& begin() {return m_range.first;}
    inline const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };

  //! square
  template<class T> inline constexpr T square( T x ) { return x*x; }

  //! radius
  template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  //! radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }

  //! calculate the average of member function called on all members in collection
  template< float (PHG4Hit::*accessor)(int) const>
  float interpolate( std::set<PHG4Hit*> hits, float rextrap )
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swr = 0;
    double swr2 = 0;
    double swx = 0;
    double swrx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x0 = (hit->*accessor)(0);
      const double x1 = (hit->*accessor)(1);
      if( std::isnan( x0 ) || std::isnan( x1 ) ) continue;

      const double w = hit->get_edep();
      if( w < 0 ) continue;

      valid = true;
      const double r0 = get_r( hit, 0 );
      const double r1 = get_r( hit, 1 );

      sw += w*2;
      swr += w*(r0 + r1);
      swr2 += w*(square(r0) + square(r1));
      swx += w*(x0 + x1);
      swrx += w*(r0*x0 + r1*x1);
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swrx - swr*swx);
    const auto beta = (swr2*swx - swr*swrx);
    const auto denom = (sw*swr2 - square(swr));

    return ( alpha*rextrap + beta )/denom;
  }

  //! true if a track is a primary
  inline int is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }


}

//_____________________________________________________________________
DSTReader::DSTReader(const std::string& name,
                     bool dryrun_info,
                     bool generateKey_info):
  SubsysReco( name),
  dryrun(dryrun_info),
  generateKey(generateKey_info)
{}

//_____________________________________________________________________
int DSTReader::Init(PHCompositeNode* topNode )
{

  //make a clusterContainer node and a trackmap node
  //assume writing macro removed them

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTReader::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "DSTReader::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }


  auto svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if( !svtxNode )
  {
    // create
    std::cout << "DSTReader::Init - SVTX node missing - creating" << std::endl;
    svtxNode = new PHCompositeNode( "SVTX" );
    dstNode->addNode(svtxNode);
  }


  auto trkrNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TRKR"));
  if( !trkrNode )
  {
    // create
    std::cout << "DSTReader::Init - TRKR node missing - creating" << std::endl;
    trkrNode = new PHCompositeNode( "TRKR" );
    dstNode->addNode(trkrNode);
  }


  auto trackNode = findNode::getClass<SvtxTrackMap>(svtxNode, "SvtxTrackMap");
  if (!trackNode)
  {
      auto newTrackNode = new PHIODataNode<PHObject>(new SvtxTrackMap_v2, "SvtxTrackMap", "PHObject");
    svtxNode->addNode(newTrackNode);
  }

  auto clsNode = findNode::getClass<TrkrClusterContainer>(trkrNode, "TRKR_CLUSTER");
  if (!clsNode)
  {
    auto newClsNode = new PHIODataNode<PHObject>(new TrkrClusterContainerv4, "TRKR_CLUSTER", "PHObject");
    trkrNode->addNode(newClsNode);
  }


  auto tpcNode = findNode::getClass<TrackSeedContainer>(svtxNode, "TpcTrackSeedContainer");
  if (!tpcNode)
  {
    auto newTpcNode = new PHIODataNode<PHObject>(new TrackSeedContainer_v1, "TpcTrackSeedContainer", "PHObject");
    svtxNode->addNode(newTpcNode);
  }

  auto siliconNode = findNode::getClass<TrackSeedContainer>(svtxNode, "SiliconTrackSeedContainer");
  if (!siliconNode)
  {
    auto newSiliconNode = new PHIODataNode<PHObject>(new TrackSeedContainer_v1, "SiliconTrackSeedContainer", "PHObject");
    svtxNode->addNode(newSiliconNode);
  }
  // evalNode->addNode(newNode);

  // auto newNode = new PHIODataNode<PHObject>( new DSTContainerv3, "DSTContainer","PHObject");
  // evalNode->addNode(newNode);
  /*
  iter = PHNodeIterator(evalNode);
  auto containerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DSTContainer"));
  if(containerNode){
        std::cout << "containerNode exists in Init"
                << "\n";
      }else{
        std::cout << "containerNode does not exists in Init"
                << "\n";
      }
    auto containerNode = dynamic_cast<DSTContainerv3*>(iter.findFirst("DSTContainer", "DSTContainer"));
  if(containerNode){
        std::cout << "containerNode exists in Init"
                << "\n";
      }else{
        std::cout << "containerNode does not exists in Init"
                << "\n";
      }
*/
  // DST container
  m_container = findNode::getClass<DSTContainerv3>(topNode, "DSTContainer");

  //m_container = findNode::getClass<DSTContainer>(topNode, "DSTContainer");
  if(!m_container){
                std::cout << "m_container does not exists in Init"
                << "\n";
      }else{
        std::cout << "m_container exists in Init"
                << "\n";
      }


  // TrackInfo container
  m_track_info_container = findNode::getClass<TrackInfoContainer_v1>(topNode, "TrackInfoContainer");

  //m_container = findNode::getClass<DSTContainer>(topNode, "DSTContainer");
  if(!m_track_info_container){
                std::cout << "m_track_info_container does not exists in Init"
                << "\n";
      }else{
        std::cout << "m_track_info_container exists in Init"
                << "\n";
      }    

  // svtxtrackmap constructer is protected
  // auto svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  // auto newNode = new PHIODataNode<PHObject>( new SvtxTrackMap, "SvtxTrackMap_2","PHObject");
  // svtxNode->addNode(newNode);

  std::cout << "init DST reading" << "\n";


  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTReader::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTReader::process_event(PHCompositeNode* topNode)
{
  std::cout << "DST event" << "\n";

  //Init(topNode);

  // load nodes
  auto res =  load_nodes(topNode);
  //std::cout <<"This is before I fill cluster map" << std::endl;
  //m_cluster_map->identify();
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  
  std::cout << "This is before evaluate track info" << std::endl;
  evaluate_track_info();
  std::cout << "This is after evaluate track info" << std::endl;
  m_track_info_container->Reset();
  


  // cleanup output container
  // if( m_container ) m_container->Reset();

  // if(m_flags&WriteEvent) 
  //   evaluate_event();
  // if(m_flags&WriteClusters) 
    //evaluate_clusters();
  //commented out evaluate clusters because it wasn't working
  //std::cout << "evaluate cluster finished" << "\n";
// if(m_flags&WriteTracks)
    //evaluate_tracks();
    //read_clusters();

    //evaluate_track_and_clusters();
  //std::cout << "evaluate track finished" << "\n";

  //std::cout <<"This is after I fill cluster map" << std::endl;
  //m_cluster_map->identify();
  // clear maps
  // m_g4hit_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTReader::End(PHCompositeNode* )
{
  std::cout << "DST reader finishes" << "\n";

  return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTReader::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_cluster_map_arr = findNode::getClass<TrkrClusterContainerv5>(topNode, "TRKR_CLUSTERV5");

  m_tpc_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  m_silicon_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<DSTContainerv3>(topNode, "DSTContainer");
  //m_container = findNode::getClass<DSTContainer>(topNode, "DSTContainer");

  m_track_info_container = findNode::getClass<TrackInfoContainer_v1>(topNode, "TrackInfoContainer");

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if(m_track_map){
        std::cout << "m_track_map exists in load"
                << "\n";
      }
      if(m_cluster_map){
        std::cout << "m_cluster_map exists in load"
                << "\n";
      }
      if(m_container){
        std::cout << "m_container exists in load"
                << "\n";
      }



  return Fun4AllReturnCodes::EVENT_OK;

}

void DSTReader::evaluate_track_and_clusters()
{

      if(m_track_map){
        std::cout << "m_track_map exists"
                << "\n";
      }
      if(m_cluster_map){
        std::cout << "m_cluster_map exists"
                << "\n";
      }
      if(m_container){
        std::cout << "m_container exists"
                << "\n";
      }

      if (!(m_track_map&&m_cluster_map&&m_container)) return;




        //there should only be as many trackContainers as layers
        std:: cout << "ArrayKeysSize: " << m_container->getArrayKeysSize() << std::endl;
        for(int iTrk = 0; iTrk < (int)m_container->getArrayKeysSize(); iTrk++){
          //auto track = (SvtxTrack_v4*)m_container->arrTrkArr->ConstructedAt(iTrk);
          unsigned int key = m_container->getArrayKeys(iTrk);
         
          //m_track_map->insert(track);
          std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << std::endl;
          auto trackContainer = (SvtxTrackArray_v1*)m_container->arrTrkContainer->ConstructedAt(iTrk);
          
          //now put this information back into original container

          //make a track and store information into it
          SvtxTrack_v4* track = new SvtxTrack_v4();
          std::cout << "track id: " << trackContainer->get_id() << std::endl;
          track->set_id(trackContainer->get_id());
          track->set_vertex_id(trackContainer->get_vertex_id());
          track->set_positive_charge(trackContainer->get_positive_charge());
          track->set_chisq(trackContainer->get_chisq());
          track->set_ndf(trackContainer->get_ndf());
          track->set_crossing(trackContainer->get_crossing());

          //insert PCA point, should be included in constructor
          //track->insert_state(new SvtxTrackState_v1(0));

          track->set_x(trackContainer->get_x());
          track->set_y(trackContainer->get_y());
          track->set_z(trackContainer->get_z());
          track->set_px(trackContainer->get_px());
          track->set_py(trackContainer->get_py());
          track->set_pz(trackContainer->get_pz());

          /*
          for( auto iter = trackContainer->begin_states(); iter != trackContainer->end_states(); ++iter )
          { track->insert_state(static_cast<SvtxTrackState*>(iter->second->CloneMe() ) ) ; }  
          */
          //Make TPCSeed and SiliconSeed
          //TrackSeed_v1* TPCSeed = new TrackSeed_v1();
          //TrackSeed_v1* SiliconSeed = new TrackSeed_v1();
          TPCSeed = new TrackSeed_v1();
          SiliconSeed = new TrackSeed_v1();

          std::cout << "Before TPC seeds" << std::endl;
          //set TPC Seed
          if(trackContainer->get_does_tpc_seed_exist()){
            TPCSeed->set_qOverR(trackContainer->tpc_seed_get_qOverR());
            TPCSeed->set_X0(trackContainer->tpc_seed_get_X0());
            TPCSeed->set_Y0(trackContainer->tpc_seed_get_Y0());
            TPCSeed->set_slope(trackContainer->tpc_seed_get_slope());
            TPCSeed->set_Z0(trackContainer->tpc_seed_get_Z0());
            TPCSeed->set_crossing(trackContainer->tpc_seed_get_crossing());

           

      
            /*
            for(int TPClayer = 7; TPClayer < 59; TPClayer++){
                if(trackContainer->getClusKey(TPClayer)==0){
              //std::cout << "No cluster at layer: " << layer << std::endl;
              continue;
            }
                
                TPCSeed->insert_cluster_key(trackContainer->getClusKey(TPClayer));
            }
            //set TPCSeed
            track->set_tpc_seed(TPCSeed);
            m_tpc_seed_container->insert(TPCSeed);
            */
          }
          std::cout << "Before silicon seeds" << std::endl;
          //set Silicon Seed
          if(trackContainer->get_does_silicon_seed_exist()){
            std::cout << "Inside silicon seeds" << std::endl;
            SiliconSeed->set_qOverR(trackContainer->silicon_seed_get_qOverR());
            SiliconSeed->set_X0(trackContainer->silicon_seed_get_X0());
            SiliconSeed->set_Y0(trackContainer->silicon_seed_get_Y0());
            SiliconSeed->set_slope(trackContainer->silicon_seed_get_slope());
            SiliconSeed->set_Z0(trackContainer->silicon_seed_get_Z0());
            SiliconSeed->set_crossing(trackContainer->silicon_seed_get_crossing());

            
            /*
            std::cout << "Before silicon seeds' clusters" << std::endl;
            for(int Siliconlayer = 0; Siliconlayer < 7; Siliconlayer++){
               
               if(trackContainer->getClusKey(Siliconlayer)==0){
              //std::cout << "No cluster at layer: " << layer << std::endl;
              continue;
            }
               SiliconSeed->insert_cluster_key(trackContainer->getClusKey(Siliconlayer));
            }
            track->set_silicon_seed(SiliconSeed);
            m_silicon_seed_container->insert(SiliconSeed);
            */
          }
          
          //add the seeds to track object
          //std::cout << "Before setting the seeds \n";
          
          
          //std::cout << "After setting the seeds \n";
          //make a cluster key set to load into this

          //delete pointers
          //delete TPCSeed;
          //delete SiliconSeed;
          //generate new clusterkey WIP

          //load it back into m_track_map
          
          //m_track_map->insert(track);


          for(int layer = 0; layer < 59; layer++){
            //check if clusterkey is default and if so then skip
            // this is 9999999999999999 in hexidecimal
            //11068046444225730969
            if(!trackContainer->getValid(layer)){
              std::cout << "No cluster at layer: " << layer << std::endl;
              continue;
            }
            //std::cout << "layer "<< layer << " clusterKey: " << clusterkey << std::endl;
            std::cout << "layer "<< layer << std::endl;
            //trackContainer->getClusKey(layer) 
            //TrkrDefs::cluskey clusterkey = trackContainer->getClusKey(layer);
            //set clusterkey
            uint8_t TrkrId = 0;
            if(layer < 3){
              TrkrId = 0;  //mvtx
            }else if(layer >= 3 && layer <7){
              TrkrId = 1; //intt
            }else{
              TrkrId = 2; //tpc
            }
            uint8_t layerbyte = (uint8_t) layer;

            std::cout << "Trkr Id: " << unsigned(TrkrId) << std::endl;
            std::cout << "Layer byte: " << unsigned(layerbyte) << std::endl;
            std::cout << "Sector Id: " << unsigned(trackContainer->getSectorId(layer)) << std::endl;
            std::cout << "Side: " << unsigned(trackContainer->getSide(layer)) << std::endl;
            

            std::cout << "Generated hitsetkey: " << ((uint32_t) TrkrId << 24) + ((uint32_t) layerbyte << 16 ) + ((uint32_t) trackContainer->getSectorId(layer) << 8) + trackContainer->getSide(layer) << std::endl;
            //(TrkrId <<24)+(layerbyte << 16 )+subsurfkey;
            //TrkrDefs::hitsetkey hitsetkey = (TrkrId << 24)+(layerbyte << 16 )+trackContainer->getSubSurfKey(layer);
            TrkrDefs::hitsetkey hitsetkey = ((uint32_t) TrkrId << 24) + ((uint32_t) layerbyte << 16 ) + ((uint32_t) trackContainer->getSectorId(layer) << 8) + trackContainer->getSide(layer);
          //make cluster object
          //TrkrClusterv5* cluster = new TrkrClusterv5();
          m_cluster = new TrkrClusterv5();

          //layer = TrkrDefs::getLayer(cluster_key);
          m_cluster->setPosition(0,trackContainer->getLocalX(layer));
          m_cluster->setPosition(1,trackContainer->getLocalY(layer));
          m_cluster->setSubSurfKey(trackContainer->getSubSurfKey(layer));
          m_cluster->setPhiError(trackContainer->getRPhiError(layer));
          m_cluster->setZError(trackContainer->getZError(layer));
          m_cluster->setAdc(trackContainer->getAdc(layer));
          m_cluster->setMaxAdc(trackContainer->getMaxAdc(layer));
          m_cluster->setPhiSize(trackContainer->getPhiSize(layer));
          m_cluster->setZSize(trackContainer->getZSize(layer));
          m_cluster->setOverlap(trackContainer->getOverlap(layer));
          m_cluster->setEdge(trackContainer->getEdge(layer));

          //TrkrDefs::hitsetkey = trackContainer->setSubSurfKey(layer);
          //load clusters back into m_cluster_map

          //key should be cluskey

          //clus key is hitsetkey for upper 32 bits
          // and cluster id for lower 32 bits
          bool foundCluster = false;
          //trackseed lists cluster keys
          for(uint32_t index = 0; index < 1000 && !foundCluster; index++){
            TrkrDefs::cluskey clusterkey = ((uint64_t)hitsetkey << 32) + index;
            
            if(!m_cluster_map->findCluster(clusterkey)){
              std::cout << "identify cluster" << std::endl;
              m_cluster->identify();
              m_cluster_map->addClusterSpecifyKey(clusterkey, m_cluster);
              foundCluster = true;



              //fill TPC and silicon seeds
              
              if(layer < 7){
                SiliconSeed->insert_cluster_key(clusterkey);
              }else{
                TPCSeed->insert_cluster_key(clusterkey);
              }
                
                
            }
            

          }
          }
          track->set_silicon_seed(SiliconSeed);
          m_silicon_seed_container->insert(SiliconSeed);
          track->set_tpc_seed(TPCSeed);
          m_tpc_seed_container->insert(TPCSeed);
          m_track_map->insertWithKey(track,key);
            //delete cluster;
          }
          //delete TPCSeed;
          //delete SiliconSeed;
        }
 




void DSTReader::evaluate_track_info(){

  for(int iTrk = 0; iTrk < (int)m_track_info_container->size(); iTrk++){ 
    std::cout << "size: " << m_track_info_container->size() << std::endl;
    std::cout << "hitbitmap: " << (m_track_info_container->get_trackinfo(iTrk))->get_hitbitmap() << std::endl;
    std::cout << "crossing: " << (m_track_info_container->get_trackinfo(iTrk))->get_crossing() << std::endl;
    std::cout << "chi^2: " << (m_track_info_container->get_trackinfo(iTrk))->get_chisq() << std::endl;
    std::cout << "NDF: " << unsigned((m_track_info_container->get_trackinfo(iTrk))->get_ndf()) << std::endl;

  }

}


//_____________________________________________________________________
void DSTReader::evaluate_event()
{
  if(!m_container) return;

  // create event struct

  // store

  // m_container->addEvent(event);
}



void DSTReader::read_clusters()
{
  std::cout << "entering read_clusters" << std::endl;
  if(!(m_cluster_map && m_cluster_map_arr&&m_hitsetcontainer)) return;

  m_cluster_map->Reset();

  std::cout << "DSTWriter::evaluate_clusters  current clusters: " << m_cluster_map_arr->size() << std::endl;
  //Int_t iCluster = 0;
  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    std::cout << "looping over hitsetkey " << std::hex << (unsigned) hitsetkey << "\n";
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map_arr->getClusters(hitsetkey)))
    {
      std::cout << "adding cluster to container" << std::endl;
      //std::cout << key << std::endl;
      cluster->identify();

      m_cluster_map->addClusterSpecifyKey(key, cluster);
      //++iCluster;
    }
    
  }

  //make loop over clustercontainerv5 and put it back into v4
  m_cluster_map_arr->Reset();
}



//_____________________________________________________________________
void DSTReader::evaluate_clusters()
{

  std::cout << "entering evaluate_clusters" << std::endl;
  // if (m_cluster_map) std::cout << "cluster map" << "\n";
  // if (m_hitsetcontainer) std::cout << "hitsetc" << "\n";
  // if (m_container) std::cout << "container" << "\n";

  // if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
  if(!(m_cluster_map && m_container)) return;

  // clear array
  // m_container->clearClusters();

  std::cout << "DSTReader::evaluate_clusters - size:" << m_cluster_map->size() << "\n";

  // std::cout << "DSTReader::evaluate_clusters - clusters: " << m_container->clusters().size() << std::endl;
  int nClusters = m_container->arrClsDST->GetEntries();
  // int nClusters = m_container->m_cluster_map_arr->size();
  std::cout << "DSTReader::evaluate_clusters - clusters (arr): " << nClusters << std::endl;
  // debugging
  // showAll();
  // showMe();
  showHitSet();
  if (!dryrun) {
    // m_cluster_map->Reset();
    std::cout << "DSTReader::evaluate_clusters - cleared. current clusters: " << m_cluster_map->size() << std::endl;
    // showMe();
    // showHitSet();
  }

  uint32_t id = 0;
  TrkrDefs::hitsetkey currentKey = 0;

  for (const auto& [hitsetkey, hitset] : range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    // auto originalRange = m_cluster_map->getClusters(hitsetkey);
    // auto iter = originalRange.first;
    // std::cout << "original vector size:" << std::distance(originalRange.second, originalRange.first) << "\n";

    for (const auto& [key, cluster] : range_adaptor(m_cluster_map_arr->getClusters(hitsetkey)))
    {
      // cluster->identify();
      std::cout << "key: " << key << std::endl;
      if (hitsetkey == currentKey)
      {
        ++id;
      }
      else
      {
        id = 0;
        currentKey = hitsetkey;
      }

      // generate cluster key
      // const auto ckey = TrkrDefs::genClusKey(key, id);
      // const auto ckeyarr = clusStruct->clusterKey;
      // std::cout << "id: " << id << ", ckey: " << ckey << ", ckeyarr: " << ckeyarr << "\n";
      // std::cout << "id: " << id <<  ", current key: " << currentKey << "\n";

      // if (ckey != ckeyarr)
      // {
      //   std::cout << "[debug] key differs!"
      //             << "\n";
      // }
      // compare with the original cluster

      auto originalRange = m_cluster_map->getClusters(hitsetkey);
      auto iter = originalRange.first;
      std::advance(iter, id);
      // std::cout << "local x of cluster: " << cluster->getLocalX() << std::endl;
      // std::cout << "local x of original: " << iter->second->getLocalX() << std::endl;
      std::cout << "comparison" << std::endl;
      if (cluster->getLocalX() != iter->second->getLocalX()) {
        std::cout << "cluster id " << id << " differs!" << std::endl;
        break;
      }
      // std::cout << "same" << std::endl;

      if (!dryrun) {
        // m_cluster_map->addClusterSpecifyKey(ckey, cluster);
        m_cluster_map->removeCluster(key);
        // std::cout << "removed key" << std::endl;
        m_cluster_map->addClusterSpecifyKey(key, cluster);
        // std::cout << "added key" << std::endl;
      }
    }
    // if (hitsetkey > 1712) {
      break;
    // }
  }
  std::cout << "DSTReader::evaluate_clusters - size:" << m_cluster_map->size() << "\n";
  // for (const auto& [hitsetkey, hitset] : range_adaptor(m_hitsetcontainer->getHitSets()))
  // {
  //   std::cout << "looping over hitsetkey " << std::hex << (unsigned) hitsetkey << std::dec << std::endl;
  //   for (const auto& [key, cluster] : range_adaptor(m_cluster_map_arr->getClusters(hitsetkey)))
  //   {
  //     cluster->identify();
  //     std::cout << key << std::endl;
  //     if (key == currentKey)
  //     {
  //       ++id;
  //     }
  //     else
  //     {
  //       id = 0;
  //       currentKey = key;
  //     }

  //     // generate cluster key
  //     // const auto ckey = TrkrDefs::genClusKey(key, id);
  //     // const auto ckeyarr = clusStruct->clusterKey;
  //     // std::cout << "id: " << id << ", ckey: " << ckey << ", ckeyarr: " << ckeyarr << "\n";
  //     // if (ckey != ckeyarr)
  //     // {
  //     //   std::cout << "[debug] key differs!"
  //     //             << "\n";
  //     // }

  //     if (!dryrun) {
  //       // m_cluster_map->addClusterSpecifyKey(ckey, cluster);
  //       // m_cluster_map->removeCluster(key);
  //       m_cluster_map->addClusterSpecifyKey(key, cluster);
  //     }
  //   }
  // }
  // std::cout << "DSTReader::evaluate_clusters - size:" << m_cluster_map->size() << "\n";

      // for (auto i = 0; i < nClusters; ++i)
      // {
      //   // DSTContainerv3::ClusterStruct clusStruct =
      //   //   *((DSTContainerv3::ClusterStruct*) m_container->arrClsDST->At(i));
      //   DSTContainerv3::ClusterStruct* clusStruct =
      //     (DSTContainerv3::ClusterStruct*) m_container->arrClsDST->At(i);
      //   // DSTContainerv3::ClusterKeyStruct clusKeyStruct =
      //   //   *((DSTContainerv3::ClusterKeyStruct*) m_container->arrKeyDST->At(i));
      //   // TrkrDefs::hitsetkey key = clusKeyStruct.hitsetkey;
      //   // TrkrClusterv4* cluster = (TrkrClusterv4*) m_container->trkrClsDST->At(i);
      //   // TrkrClusterv4* cluster = (TrkrClusterv4*) m_container->trkrClsDST->At(i);
      //   // std::cout << "cluster " << (unsigned) id << "\n";
      //   // TrkrCluster newCluster = recover_cluster(cluster);

      //   // auto newCluster = clusStruct.cluster;
      //   // auto key = clusStruct.hitsetkey;
      //   // auto newCluster = &clusStruct->cluster;
      //   auto newCluster = new TrkrClusterv4();
      //   newCluster->CopyFrom(&clusStruct->cluster);

      //   auto key = clusStruct->hitsetkey;
      //   // auto newCluster = new TrkrClusterv4;
      //   // newCluster->setLocalX(clusStruct.loc_x);
      //   // newCluster->setLocalY(clusStruct.loc_y);
      //   // // TrkrDefs::hitsetkey key = 0;
      //   // // // std::cout << std::bitset<32>(key) << "\n";
      //   // // // unsigned char z_unsigned = (unsigned) clusStruct.z_seg;
      //   // // // key |= (z_unsigned << TrkrDefs::kBitShiftZElement);
      //   // // key |= (clusStruct.z_seg << TrkrDefs::kBitShiftZElement);
      //   // // // std::cout << std::bitset<32>(key) << "\n";
      //   // // // unsigned char phi_unsigned = (unsigned) clusStruct.phi_seg;
      //   // // // key |= (phi_unsigned << TrkrDefs::kBitShiftPhiElement);
      //   // // key |= (clusStruct.phi_seg << TrkrDefs::kBitShiftPhiElement);
      //   // // // std::cout << std::bitset<32>(key) << "\n";
      //   // // // unsigned char layer_unsigned = (unsigned) clusStruct.layer;
      //   // // // key |= (layer_unsigned << TrkrDefs::kBitShiftLayer);
      //   // // key |= (clusStruct.layer << TrkrDefs::kBitShiftLayer);
      //   // // // std::cout << std::bitset<32>(key) << "\n";
      //   // // key |= (id << TrkrDefs::kBitShiftTrkrId);
      //   // // // std::cout << std::bitset<32>(key) << "\n";
      //   // // TrkrDefs::cluskey cluskey = key;
      //   // // cluskey = (cluskey << TrkrDefs::kBitShiftClusId);
      //   // // // std::cout << std::bitset<64>(cluskey) << "\n";

      //   // // set size and error
      //   // // for (int j = 0; j < 3; ++j) {
      //   // //   for (int i = 0; i < 3; ++i) {
      //   // //     newCluster->setSize(i, j, DSTContainerv3::covarIndex(i, j));
      //   // //     newCluster->setError(i, j, DSTContainerv3::covarIndex(i, j));
      //   // //   }
      //   // // }

      //   // int nLocal = 2;
      //   // for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      //   //   for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
      //   //     newCluster->setActsLocalError(iLocal, jLocal,
      //   //                                  clusStruct.actsLocalError[iLocal][jLocal]);
      //   //     // std::cout << "actslocalerror:" << newCluster->
      //   //       // getActsLocalError(iLocal, jLocal) << "\n";
      //   //   }
      //   // }
      //   // newCluster->setSubSurfKey(clusStruct.subSurfKey);
      //   // // std::cout << "subsurfkey: " << newCluster->getSubSurfKey() << "\n";

      //   // newCluster->setAdc(clusStruct.adc);
      //   // // std::cout << "adc: " << newCluster->getAdc() << "\n";
      //   // m_cluster_map->addClusterSpecifyKey(clusStruct.clusterKey, newCluster);

      //   newCluster->identify();
      //   // newCluster.identify();
      //   std::cout << key << "\n";
      //   if (key == currentKey) {
      //     ++id;
      //   } else {
      //     id = 0;
      //     currentKey = key;
      //   }

      //   // std::cout << "print" << "\n";
      //   // printCluster(*newCluster);

      //   // cluster->identify();
      //   // generate cluster key
      //   const auto ckey = TrkrDefs::genClusKey( key, id );
      //   const auto ckeyarr = clusStruct->clusterKey;
      //   std::cout << "id: " << id << ", ckey: " << ckey << ", ckeyarr: " << ckeyarr << "\n";
      //   if (ckey != ckeyarr) {
      //     std::cout << "[debug] key differs!" << "\n";
      //   }

      //   // m_cluster_map->addClusterSpecifyKey(clusStruct.clusterKey, newCluster);
      //   // m_cluster_map->addClusterSpecifyKey(ckey, cluster);
      //   if (!dryrun) {
      //     if (generateKey) {
      //       m_cluster_map->removeCluster(ckey);
      //       m_cluster_map->addClusterSpecifyKey(ckey, newCluster);
      //     } else {
      //       m_cluster_map->removeCluster(ckeyarr);
      //       m_cluster_map->addClusterSpecifyKey(ckeyarr, newCluster);
      //     }
      //   }
      //   // ++id;
      // }

      // showMe();
      // showHitSet();
      // for (auto& clusStruct : m_container->clusters())
      // {
      //   std::cout << "cluster " << (unsigned) id << "\n";
      //   // TrkrCluster newCluster = recover_cluster(cluster);

      //   auto newCluster = new TrkrClusterv3;
      //   newCluster->setLocalX(clusterStruct.loc_x);
      //   newCluster->setLocalY(clusterStruct.loc_y);
      //   TrkrDefs::hitsetkey key = 0;
      //   // std::cout << std::bitset<32>(key) << "\n";
      //   // unsigned char z_unsigned = (unsigned) clusterStruct.z_seg;
      //   // key |= (z_unsigned << TrkrDefs::kBitShiftZElement);
      //   key |= (clusterStruct.z_seg << TrkrDefs::kBitShiftZElement);
      //   // std::cout << std::bitset<32>(key) << "\n";
      //   // unsigned char phi_unsigned = (unsigned) clusterStruct.phi_seg;
      //   // key |= (phi_unsigned << TrkrDefs::kBitShiftPhiElement);
      //   key |= (clusterStruct.phi_seg << TrkrDefs::kBitShiftPhiElement);
      //   // std::cout << std::bitset<32>(key) << "\n";
      //   // unsigned char layer_unsigned = (unsigned) clusterStruct.layer;
      //   // key |= (layer_unsigned << TrkrDefs::kBitShiftLayer);
      //   key |= (clusterStruct.layer << TrkrDefs::kBitShiftLayer);
      //   // std::cout << std::bitset<32>(key) << "\n";
      //   key |= (id << TrkrDefs::kBitShiftTrkrId);
      //   // std::cout << std::bitset<32>(key) << "\n";
      //   TrkrDefs::cluskey cluskey = key;
      //   cluskey = (cluskey << TrkrDefs::kBitShiftClusId);
      //   // std::cout << std::bitset<64>(cluskey) << "\n";

      //   // set size and error
      //   // for (int j = 0; j < 3; ++j) {
      //   //   for (int i = 0; i < 3; ++i) {
      //   //     newCluster->setSize(i, j, DSTContainerv3::covarIndex(i, j));
      //   //     newCluster->setError(i, j, DSTContainerv3::covarIndex(i, j));
      //   //   }
      //   // }

      //   int nLocal = 2;
      //   for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      //     for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
      //       newCluster->setActsLocalError(iLocal, jLocal,
      //                                    clusterStruct.actsLocalError[iLocal][jLocal]);
      //       // std::cout << "actslocalerror:" << newCluster->
      //         // getActsLocalError(iLocal, jLocal) << "\n";
      //     }
      //   }
      //   newCluster->setSubSurfKey(clusterStruct.subSurfKey);
      //   // std::cout << "subsurfkey: " << newCluster->getSubSurfKey() << "\n";

      //   newCluster->setAdc(clusterStruct.adc);
      //   // std::cout << "adc: " << newCluster->getAdc() << "\n";

      //   // debugging
      //   // std::cout << "key from clust: " << std::hex << clusterStruct.clusterKey << "\n";
      //   // std::cout << "key by me: " << std::hex << cluskey << std::dec << "\n";
      //   // for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      //   //   std::cout << "positions: " << newCluster->getPosition(iLocal) << "\n";
      //   // }
      //   // std::cout << "rphierror: " << newCluster->getRPhiError() << "\n";
      //   // std::cout << "zerror: " << newCluster->getZError() << "\n";

      //   // get hitsetkey from cluster
      //   // const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( cluskey );
      //   // // get cluster index in vector
      //   // const auto index = TrkrDefs::getClusIndex( cluskey );
      //   // std::cout << "index:" << index << "\n";

      //   // m_cluster_map->addClusterSpecifyKey(cluskey, &newCluster);
      //   m_cluster_map->addClusterSpecifyKey(clusterStruct.clusterKey, newCluster);
      //   ++id;
      //   std::cout << "DSTReader::evaluate_clusters - saved clusters: " << m_cluster_map->size() << std::endl;
      //   showMe();
      // }

      std::cout << "DSTReader::evaluate_clusters - saved clusters: " << m_cluster_map->size() << std::endl;

      // showAll();
    }

    //_____________________________________________________________________
    void DSTReader::evaluate_tracks()
    {
      if(m_track_map){
        std::cout << "m_track_map exists"
                << "\n";
      }
      if(m_cluster_map){
        std::cout << "m_cluster_map exists"
                << "\n";
      }
      if(m_container){
        std::cout << "m_container exists"
                << "\n";
      }

      if (!(m_track_map && m_cluster_map && m_container)) return;

      std::cout << "reading file"
                << "\n";

      std::cout << "DSTReader: " << m_container->tracks().size() << std::endl;

      // // clear array
      // m_container->clearTracks();

      // for( const auto& trackpair:*m_track_map )
      // {

      //   const auto track = trackpair.second;
      //   auto track_struct = create_track( track );

      //   // truth information
      //   const auto [id,contributors] = get_max_contributor( track );
      //   track_struct.contributors = contributors;

      //   // get particle
      //   auto particle = m_g4truthinfo->GetParticle(id);
      //   track_struct.embed = get_embed(particle);
      //   add_truth_information(track_struct, particle);

      //   // running iterator over track states, used to match a given cluster to a track state
      //   auto state_iter = track->begin_states();

      //   // loop over clusters
      //   for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
      //   {

      //     const auto& cluster_key = *key_iter;
      //     auto cluster = m_cluster_map->findCluster( cluster_key );
      //     if( !cluster )
      //     {
      //       std::cout << "DSTReader::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
      //       continue;
      //     }

      //     // create new cluster struct
      //     auto cluster_struct = create_cluster( cluster_key, cluster );
      //     add_cluster_size( cluster_struct, cluster_key, m_cluster_hit_map );
      //     add_cluster_energy( cluster_struct, cluster_key, m_cluster_hit_map, m_hitsetcontainer );

      //     // truth information
      //     const auto g4hits = find_g4hits( cluster_key );
      //     add_truth_information( cluster_struct, g4hits );

      //     // find track state that is the closest to cluster
      //     // this assumes that both clusters and states are sorted along r
      //     const auto radius( cluster_struct.r );
      //     float dr_min = -1;
      //     for( auto iter = state_iter; iter != track->end_states(); ++iter )
      //     {
      //       const auto dr = std::abs( radius - get_r( iter->second->get_x(), iter->second->get_y() ) );
      //       if( dr_min < 0 || dr < dr_min )
      //       {
      //         state_iter = iter;
      //         dr_min = dr;
      //       } else break;
      //     }

      //     // store track state in cluster struct
      //     add_trk_information( cluster_struct, state_iter->second );

      //     // add to track
      //     track_struct.clusters.push_back( cluster_struct );
      //   }
      // m_container->addTrack( track_struct );
      // }

      std::cout << "DSTReader::evaluate_tracks - tracks: " << m_container->tracks().size() << std::endl;
      m_track_map->Reset();

    /*
      for (auto trackStruct : m_container->tracks())
      {
        SvtxTrack_v1 track;
        track.set_charge(trackStruct.charge);

        track.set_chisq(trackStruct.chisquare);
        track.set_ndf(trackStruct.ndf);

        track.set_x(trackStruct.x);
        track.set_y(trackStruct.y);
        track.set_z(trackStruct.z);

        track.set_px(trackStruct.px);
        track.set_py(trackStruct.py);
        track.set_pz(trackStruct.pz);

        m_track_map->insert(&track);
      }
*/
      //Set a track from TrackArray similar to how this sets a track for TrackMap
        //SvtxTrack_v4 track;
        //m_container->arrTrkArr;
        /*
        for(int iTrk = 0; iTrk < m_container->arrTrkArr->Size(); iTrk++){
          auto track = m_container->arrTrkArr->ConstructedAt(iTrk);
          unsigned int key = m_container->TrackArrayKeys[iTrk];
          m_track_map->insertWithKey(&track,key);
        }
        */
        std:: cout << "ArrayKeysSize: " << m_container->getArrayKeysSize() << std::endl;
        for(int iTrk = 0; iTrk < (int)m_container->getArrayKeysSize(); iTrk++){
          auto track = (SvtxTrack_v4*)m_container->arrTrkArr->ConstructedAt(iTrk);
          unsigned int key = m_container->getArrayKeys(iTrk);
          m_track_map->insertWithKey(track,key);
          //m_track_map->insert(track);
        }

      std::cout << "DSTReader::evaluate_tracks - saved tracks: " << m_track_map->size() << std::endl;


      


    }

    //_____________________________________________________________________
    DSTReader::G4HitSet DSTReader::find_g4hits(TrkrDefs::cluskey cluster_key) const
    {
      // check maps
      if (!(m_cluster_hit_map && m_hit_truth_map)) return G4HitSet();

      // check if in map
      auto map_iter = m_g4hit_map.lower_bound(cluster_key);
      if (map_iter != m_g4hit_map.end() && cluster_key == map_iter->first)
      {
        return map_iter->second;
      }

      // find hitset associated to cluster
      G4HitSet out;
      const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
      for (const auto& [first, hit_key] : range_adaptor(m_cluster_hit_map->getHits(cluster_key)))
      {
        // store hits to g4hit associations
        TrkrHitTruthAssoc::MMap g4hit_map;
        m_hit_truth_map->getG4Hits(hitset_key, hit_key, g4hit_map);

        // find corresponding g4 hist
        for (auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter)
        {
          const auto g4hit_key = truth_iter->second.second;
          PHG4Hit* g4hit = nullptr;

          switch (TrkrDefs::getTrkrId(hitset_key))
          {
          case TrkrDefs::mvtxId:
            if (m_g4hits_mvtx) g4hit = m_g4hits_mvtx->findHit(g4hit_key);
            break;

          case TrkrDefs::inttId:
            if (m_g4hits_intt) g4hit = m_g4hits_intt->findHit(g4hit_key);
            break;

          case TrkrDefs::tpcId:
            if (m_g4hits_tpc) g4hit = m_g4hits_tpc->findHit(g4hit_key);
            break;

          case TrkrDefs::micromegasId:
            if (m_g4hits_micromegas) g4hit = m_g4hits_micromegas->findHit(g4hit_key);
            break;

          default:
            break;
          }

          if (g4hit)
            out.insert(g4hit);
          else
            std::cout << "DSTReader::find_g4hits - g4hit not found " << g4hit_key << std::endl;
        }
      }

      // insert in map and return
      return m_g4hit_map.insert(map_iter, std::make_pair(cluster_key, std::move(out)))->second;
    }

    //_____________________________________________________________________
    std::pair<int, int> DSTReader::get_max_contributor(SvtxTrack * track) const
    {
      if (!(m_track_map && m_cluster_map && m_g4truthinfo)) return {0, 0};

      // maps MC track id and number of matching g4hits
      using IdMap = std::map<int, int>;
      IdMap contributor_map;

      // loop over clusters
      for (auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        for (const auto& hit : find_g4hits(cluster_key))
        {
          const int trkid = hit->get_trkid();
          auto iter = contributor_map.lower_bound(trkid);
          if (iter == contributor_map.end() || iter->first != trkid)
          {
            contributor_map.insert(iter, std::make_pair(trkid, 1));
          }
          else
            ++iter->second;
        }
      }

      if (contributor_map.empty())
        return {0, 0};
      else
        return *std::max_element(
            contributor_map.cbegin(), contributor_map.cend(),
            [](const IdMap::value_type& first, const IdMap::value_type& second)
            { return first.second < second.second; });
    }

    // //! create svx track from struct
    // SvtxTrack DSTReader::recover_track(DSTContainerv3::TrackStruct trackStruct)
    // {
    //   SvtxTrack_v1 track;
    //   track.set_charge(trackStruct.charge);
    //   // trackStruct.nclusters = track->size_cluster_keys();
    //   // trackStruct.mask = get_mask( track );

    //   // track.insert_cluster_key();
    //   // trackStruct.nclusters_mvtx = get_clusters<TrkrDefs::mvtxId>( track );
    //   // trackStruct.nclusters_intt = get_clusters<TrkrDefs::inttId>( track );
    //   // trackStruct.nclusters_tpc = get_clusters<TrkrDefs::tpcId>( track );
    //   // trackStruct.nclusters_micromegas = get_clusters<TrkrDefs::micromegasId>( track );

    //   track.set_chisq(trackStruct.chisquare);
    //   track.set_ndf(trackStruct.ndf);

    //   track.set_x(trackStruct.x);
    //   track.set_y(trackStruct.y);
    //   track.set_z(trackStruct.z);

    //   track.set_px(trackStruct.px);
    //   track.set_py(trackStruct.py);
    //   track.set_pz(trackStruct.pz);

    //   return track;
    // }

    // TrkrCluster DSTReader::recover_cluster(DSTContainerv3::ClusterStruct clusterStruct)
    // {
    //   TrkrClusterv3 cluster;
    //   cluster.setLocalX(clusterStruct.loc_x);
    //   cluster.setLocalY(clusterStruct.loc_y);

    //   return cluster;
    // }

    //_____________________________________________________________________
    int DSTReader::get_embed(PHG4Particle * particle) const
    {
      return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded(particle->get_primary_id()) : 0;
    }

    void DSTReader::printCluster(TrkrCluster & newCluster) const
    {
      int nLocal = 2;
      for (auto iLocal = 0; iLocal < nLocal; ++iLocal)
      {
        for (auto jLocal = 0; jLocal < nLocal; ++jLocal)
        {
          std::cout << "actslocalerror:" << newCluster.getActsLocalError(iLocal, jLocal) << std::endl;
        }
      }
      std::cout << "subsurfkey: " << newCluster.getSubSurfKey() << std::endl;
      std::cout << "adc: " << newCluster.getAdc() << std::endl;
      for (auto iLocal = 0; iLocal < nLocal; ++iLocal)
      {
        std::cout << "positions: " << newCluster.getPosition(iLocal) << std::endl;
      }
      std::cout << "rphierror: " << newCluster.getRPhiError() << std::endl;
      std::cout << "zerror: " << newCluster.getZError() << std::endl;
    }

    void DSTReader::showMe() const
    {
      std::cout << "hitset size: " << m_hitsetcontainer->size() << std::endl;
      for (const auto& [hitsetkey, hitset] : range_adaptor(m_hitsetcontainer->getHitSets()))
      {
        std::cout << "looping over hitsetkey " << std::hex << hitsetkey << std::endl;
        for (const auto& [key, cluster] : range_adaptor(m_cluster_map->getClusters(hitsetkey)))
        {
          std::cout << "printCluster " << std::hex << key << std::dec << std::endl;
          // printCluster(*cluster);
          // std::cout << "cluster is valid " << cluster->isValid() << std::endl;
          cluster->identify();
          std::cout << "loc x: " << cluster->getLocalX() << ", loc y: " << cluster->getLocalY() << ", subsurfkey: " << cluster->getSubSurfKey() << ", adc: " << cluster->getAdc() << ", phisize: " << cluster->getPhiSize() << ", zsize: " << cluster->getZSize() << "\n";
        }
        // break;
      }
    }

    void DSTReader::showHitSet() const
    {
      std::cout << "hitset size: " << m_hitsetcontainer->size() << std::endl;
      for (const auto& [hitsetkey, hitset] : range_adaptor(m_hitsetcontainer->getHitSets()))
      {
        std::cout << "looping over hitsetkey " << std::hex << hitsetkey << std::dec << std::endl;
        // for (const auto& [key, cluster] : range_adaptor(m_cluster_map->getClusters(hitsetkey)))
        // {
      }
    }

    void DSTReader::showAll() const
    {
      std::cout << "show all clusters of " << m_cluster_map->size() << std::endl;
      for (const auto& hitset : m_cluster_map->getHitSetKeys())
      {
        for (const auto& [key, cluster] : range_adaptor(m_cluster_map->getClusters(hitset)))
        {
          std::cout << "identify cluster with key " << std::hex << key << std::dec << std::endl;
          std::cout << "cluster is valid? " << cluster->isValid() << std::endl;
          cluster->identify();
          std::cout << "our time"
                    << "\n";
          printCluster(*cluster);
        }
      }
    }

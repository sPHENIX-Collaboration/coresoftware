/*!
 * \file DSTWriter.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "DSTWriter.h"
#include "DSTContainerv3.h"
#include "DSTContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <trackbase/InttDefs.h>
#include <micromegas/MicromegasDefs.h>
#include <trackbase/MvtxDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
//include new cluster
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrClusterContainer.h>
// #include <trackbase/TrkrClusterContainerv5.h>
#include <trackbase/TrkrClusterContainerv5.h>
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
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackArray_v1.h>
#include <trackbase_historic/SvtxTrackInfo_v1.h>
#include <trackbase_historic/TrackInfoContainer_v1.h>
#include <trackbase_historic/TrackStateInfo_v1.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
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

  //! get mask from track clusters
  int64_t get_mask( SvtxTrack* track )
  { return std::accumulate( track->begin_cluster_keys(), track->end_cluster_keys(), int64_t(0),
      []( int64_t value, const TrkrDefs::cluskey& key ) {
        return TrkrDefs::getLayer(key)<64 ? value|(1LL<<TrkrDefs::getLayer(key)) : value;
      } );
  }

  //! return number of clusters of a given type
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    return std::count_if( track->begin_cluster_keys(), track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }

  //! create track struct from struct from svx track
  DSTContainerv3::TrackStruct create_track( SvtxTrack* track )
  {
    DSTContainerv3::TrackStruct trackStruct;

    trackStruct.charge = track->get_charge();
    trackStruct.nclusters = track->size_cluster_keys();
    trackStruct.mask = get_mask( track );
    trackStruct.nclusters_mvtx = get_clusters<TrkrDefs::mvtxId>( track );
    trackStruct.nclusters_intt = get_clusters<TrkrDefs::inttId>( track );
    trackStruct.nclusters_tpc = get_clusters<TrkrDefs::tpcId>( track );
    trackStruct.nclusters_micromegas = get_clusters<TrkrDefs::micromegasId>( track );

    trackStruct.chisquare = track->get_chisq();
    trackStruct.ndf = track->get_ndf();

    trackStruct.x = track->get_x();
    trackStruct.y = track->get_y();
    trackStruct.z = track->get_z();

    trackStruct.px = track->get_px();
    trackStruct.py = track->get_py();
    trackStruct.pz = track->get_pz();

    return trackStruct;
  }

  // //! create cluster struct from svx cluster
  // DSTContainerv3::ClusterStruct create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster )
  // {
  //   DSTContainerv3::ClusterStruct cluster_struct;
  //   cluster_struct.clusterKey = key;
  //   cluster_struct.layer = TrkrDefs::getLayer(key);
  //   cluster_struct.phi_seg = TrkrDefs::getPhiElement(key);
  //   cluster_struct.z_seg = TrkrDefs::getZElement(key);
  //   cluster_struct.loc_x = cluster->getLocalX();
  //   cluster_struct.loc_y = cluster->getLocalY();
  //   cluster_struct.loc_x_err = cluster->getActsLocalError(0,0);
  //   cluster_struct.loc_y_err = cluster->getActsLocalError(1,1);
  //   cluster_struct.adc = cluster->getAdc();
  //   // for (int j = 0; j < 3; ++j) {
  //   //   for (int i = 0; i < 3; ++i) {
  //   //     cluster_struct.cor_size[DSTContainerv3::covarIndex(i, j)] = cluster->getSize(i, j);
  //   //     cluster_struct.cor_error[DSTContainerv3::covarIndex(i, j)] = cluster->getError(i, j);
  //   //   }
  //   // }

  //   // for v3
  //   int nLocal = 2;
  //   for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
  //     for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
  //       cluster_struct.actsLocalError[iLocal][jLocal] =
  //         cluster->getActsLocalError(iLocal, jLocal);
  //     }
  //   }

  //   cluster_struct.subSurfKey = cluster->getSubSurfKey();


  //   // below for v4
  //   // cluster_struct.overlap = cluster->getOverlap();
  //   // cluster_struct.edge = cluster->getEdge();

  //   return cluster_struct;
  // }

  //! number of hits associated to cluster
  // void add_cluster_size( DSTContainerv3::ClusterStruct& cluster, TrkrDefs::cluskey clus_key, TrkrClusterHitAssoc* cluster_hit_map )
  // {
  //   if( !cluster_hit_map ) return;
  //   const auto range = cluster_hit_map->getHits(clus_key);

  //   // store full size
  //   cluster.size =  std::distance( range.first, range.second );

  //   const auto detId = TrkrDefs::getTrkrId(clus_key);
  //   if(detId == TrkrDefs::micromegasId)
  //   {

  //     // for micromegas the directional cluster size depends on segmentation type
  //     auto segmentation_type = MicromegasDefs::getSegmentationType(clus_key);
  //     if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_Z ) cluster.z_size = cluster.size;
  //     else cluster.phi_size = cluster.size;

  //   } else {

  //     // for other detectors, one must loop over the constituting hits
  //     std::set<int> phibins;
  //     std::set<int> zbins;
  //     for(const auto& [first, hit_key]:range_adaptor(range))
  //     {
  //       switch( detId )
  //       {
  //         default: break;
  //         case TrkrDefs::mvtxId:
  //         {
  //           phibins.insert( MvtxDefs::getRow( hit_key ) );
  //           zbins.insert( MvtxDefs::getCol( hit_key ) );
  //           break;
  //         }
  //         case TrkrDefs::inttId:
  //         {
  //           phibins.insert( InttDefs::getRow( hit_key ) );
  //           zbins.insert( InttDefs::getCol( hit_key ) );
  //           break;
  //         }
  //         case TrkrDefs::tpcId:
  //         {
  //           phibins.insert( TpcDefs::getPad( hit_key ) );
  //           zbins.insert( TpcDefs::getTBin( hit_key ) );
  //           break;
  //         }
  //       }
  //     }
  //     cluster.phi_size = phibins.size();
  //     cluster.z_size = zbins.size();
  //   }
  // }


  // //! hit energy for a given cluster
  // void add_cluster_energy( DSTContainerv3::ClusterStruct& cluster, TrkrDefs::cluskey clus_key,
  //   TrkrClusterHitAssoc* cluster_hit_map,
  //   TrkrHitSetContainer* hitsetcontainer )
  // {

  //   // check container
  //   if(!(cluster_hit_map && hitsetcontainer)) return;

  //   // for now this is only filled for micromegas
  //   const auto detId = TrkrDefs::getTrkrId(clus_key);
  //   if(detId != TrkrDefs::micromegasId) return;

  //   const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(clus_key);
  //   const auto hitset = hitsetcontainer->findHitSet( hitset_key );
  //   if( !hitset ) return;

  //   const auto range = cluster_hit_map->getHits(clus_key);
  //   float energy_max = 0;
  //   float energy_sum = 0;

  //   for( const auto& pair:range_adaptor(range))
  //   {
  //     const auto hit = hitset->getHit( pair.second );
  //     if( hit )
  //     {
  //       const auto energy = hit->getEnergy();
  //       energy_sum += energy;
  //       if( energy > energy_max ) energy_max = energy;
  //     }
  //   }
  //   cluster.amp = energy_sum;
  // }

}

//_____________________________________________________________________
DSTWriter::DSTWriter( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int DSTWriter::Init(PHCompositeNode* topNode )
{
  std::cout << "Writer Init start" << std::endl;
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTWriter::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "DSTWriter::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  // // TClonesArary
  // m_container->arrClsDST = new TClonesArray("DSTContainerv3::ClusterStruct");
  // m_container->trkrClsDST = new TClonesArray("TrkrClusterv4");

  auto newNode = new PHIODataNode<PHObject>( new DSTContainerv3, "DSTContainer","PHObject");
  evalNode->addNode(newNode);

  auto newInfoNode = new PHIODataNode<PHObject>( new TrackInfoContainer_v1, "TrackInfoContainer","PHObject");
  evalNode->addNode(newInfoNode);

  auto newClsNode = new PHIODataNode<PHObject>(new TrkrClusterContainerv5, "TRKR_CLUSTERV5", "PHObject");
  evalNode->addNode(newClsNode);

  //auto newTrackNode = new PHIODataNode<PHObject>(new SvtxTrackArray_v1, "TRACK_ARRAYV1", "PHObject");
  //evalNode->addNode(newTrackNode);

  // TClonesArary
  // fcl = new TFile("dstcl.root", "recreate");
  // tcl = new TTree("tcl", "dst tree");
  // arrEvt = new TClonesArray("DSTContainerv3::EventStruct");
  // arrTrk = new TClonesArray("DSTContainerv3::TrackStruct");
  // arrCls = new TClonesArray("DSTContainerv3::ClusterStruct");
  // tcl->Branch("evt", &arrEvt);
  // tcl->Branch("trk", &arrTrk);
  // tcl->Branch("cls", &arrCls);
  std::cout << "Writer Init end" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTWriter::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTWriter::process_event(PHCompositeNode* topNode)
{
  //make topNode run in Init
  //Init(topNode);
  // load nodes
  std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << std::endl;
  std::cout << "DSTWriter::process_event" << std::endl;
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;
  std::cout << "Return codes  start" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  // cleanup output container
  if( m_container ) m_container->Reset();
  if( m_track_info_container) m_track_info_container->Reset();
  


      



  //if(m_flags&WriteEvent)
    //evaluate_event();

  //comment out WriteTracks&WriteClusters for new method
  
  //if(m_flags&WriteClusters) {
    //std::cout << "DSTWriter::process_event doing clusters" << std::endl;
    //evaluate_clusters();
  //}
  /*
  if(m_flags&WriteTracks) {
    std::cout << "DSTWriter::process_event doing tracks" << std::endl;
    evaluate_tracks();
  }
  */
 /*
if(m_flags){
  std::cout << "DSTWriter::m_flags is true" << std::endl;
}
if(WriteTracks){
  std::cout << "DSTWriter::WriteTracks is true" << std::endl;
}
if(WriteClusters){
  std::cout << "DSTWriter::WriteClusters is true" << std::endl;
}
*/
  //if(m_flags&WriteTracks&WriteClusters) {
    //std::cout << "DSTWriter::process_event doing tracks" << std::endl;
    //evaluate_track_and_cluster();
  //}



  std::cout << "Evalutate track info" << std::endl;
  evaluate_track_info();

  std::cout << "exiting event" << "\n";


  // clear maps
  std::cout << "Before map clear" << std::endl;
  m_g4hit_map.clear();
  std::cout << "After map clear" << std::endl;
  std::cout << "Return codes end" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTWriter::End(PHCompositeNode* )
{
  std::cout << "In the end" << std::endl;
  // tcl->Print();
  // tcl->Write();
  // fcl->Write();
  // fcl->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTWriter::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  //added array node
  //don't need node because this will be stored in m_container
  //m_track_array = findNode:getClass<SvtxTrackArray_v1>(topNode, "TRACK_ARRAYV1");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_cluster_map_arr = findNode::getClass<TrkrClusterContainerv5>(topNode, "TRKR_CLUSTERV5");
  
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

  // svtx trackseed container
  m_svtxtrackseed = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");


  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return Fun4AllReturnCodes::EVENT_OK;

}

//____________________________________________________________________
void DSTWriter::evaluate_track_info(){


  if(!(m_track_info_container)) return;

  m_track_info_container->Reset();
  //get track into track info
  Int_t iTrk = 0;
  //long unsigned int iKey = 0;


  std::cout << "Before loop" << "\n";
  for( const auto& trackpair:*m_track_map )
  {
    const auto track = trackpair.second;
    //this track will have a TPC and Silicon seed

    uint64_t hitbitmap = 0; 

    SvtxTrackInfo_v1* trackInfo = new SvtxTrackInfo_v1();

      std::cout << "Before seeds" << "\n";
    TrackSeed* TPCSeed = track->get_tpc_seed();
    TrackSeed* SiliconSeed = track->get_silicon_seed();
    if(TPCSeed){
      for( auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter )
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_cluster_map->findCluster( cluster_key );
        if( !cluster )
        {
          std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        //store information in track array
        //std::cout << "TPC clusterkey: " << cluster_key <<"\n";
        //std::cout << "TPC subsurfkey: " << cluster->getSubSurfKey() << std::endl;
        uint8_t layer = TrkrDefs::getLayer(cluster_key);
        std::cout << "Layer is: " << unsigned(layer) << std::endl;
        hitbitmap = hitbitmap + ((uint64_t)1 << layer);
        

        //TrkrDefs::
      
      }
    }
    std::cout << "Before Silicon seeds" << "\n";

    if(!SiliconSeed){
      std::cout << "Silicon Seed does not exist" << std::endl;
    }

    if(SiliconSeed){
      for( auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter )
      {
        const auto& cluster_key = *key_iter;
       auto cluster = m_cluster_map->findCluster( cluster_key );
       if( !cluster )
       {
         std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
         continue;
        }
        uint8_t layer = TrkrDefs::getLayer(cluster_key);
        std::cout << "Layer is: " << unsigned(layer) << std::endl;
        hitbitmap = hitbitmap + ((uint64_t)1 << layer);
     }
   }
/*
    std::cout << "After Track seeds" << "\n";
    trackInfo.set_hitbitmap(hitbitmap);
    trackInfo.set_chisq(track->get_chisq());
    trackInfo.set_ndf(track->get_ndf());


    trackInfo.set_x(track->get_x());
    trackInfo.set_y(track->get_y());
    trackInfo.set_z(track->get_z());
    trackInfo.set_px(track->get_px());
    trackInfo.set_py(track->get_py());
    trackInfo.set_pz(track->get_pz());

    std::cout << "track.get_z(): " << track->get_z() << std::endl;
    std::cout << "trackInfo.get_z(): " << trackInfo.get_z() << std::endl;

    int covarianceIndex = 0;
    for(int i = 0; i<6; i++){
      for(int j = i; j<6; j++){
        std::cout << "covariance index: " << covarianceIndex << std::endl;
        trackInfo.set_covariance(covarianceIndex, track->get_error(i,j));
        covarianceIndex++;
      }
    }
*/
  std::cout << "After Track seeds" << "\n";
    trackInfo->set_hitbitmap(hitbitmap);
    trackInfo->set_chisq(track->get_chisq());
    trackInfo->set_ndf(track->get_ndf());
    trackInfo->set_crossing(track->get_crossing());

    //this is a test
    //trackInfo->set_crossing(5);  

    trackInfo->set_x(track->get_x());
    trackInfo->set_y(track->get_y());
    trackInfo->set_z(track->get_z());
    trackInfo->set_px(track->get_px());
    trackInfo->set_py(track->get_py());
    trackInfo->set_pz(track->get_pz());

    std::cout << "track crossing: " <<track->get_crossing() << std::endl;

    std::cout << "track.get_z(): " << track->get_z() << std::endl;
    std::cout << "trackInfo.get_z(): " << trackInfo->get_z() << std::endl;
    std::cout << "Hitbitmap: " << trackInfo->get_hitbitmap() << std::endl;
    std::cout << "crossing: " << trackInfo->get_crossing() << std::endl;
    std::cout << "chi^2: " << trackInfo->get_chisq() << std::endl;
    std::cout << "ndf: " << unsigned(trackInfo->get_ndf()) << std::endl;

    int covarianceIndex = 0;
    for(int i = 0; i<6; i++){
      for(int j = i; j<6; j++){
        std::cout << "covariance index: " << covarianceIndex << std::endl;
        trackInfo->set_covariance(covarianceIndex, track->get_error(i,j));
        covarianceIndex++;
      }
    }



    std::cout << "Right before adding track info" << iTrk << std::endl;
    m_track_info_container->add_trackinfo(iTrk, trackInfo);
    std::cout << "Right after adding track info" << std::endl;
    delete trackInfo;
    ++iTrk;
  }

  


  //add trackinfo to trackinfocontainer

}
//_____________________________________________________________________

//_____________________________________________________________________
void DSTWriter::evaluate_track_and_cluster(){
//use this to create object that looks through both tracks and clusters and saves into new object
//make sure clusters exist
//std::cout << "start of check" << "\n";
if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
//make sure tracks exist
if( !( m_track_map && m_cluster_map && m_container ) ) return;
//std::cout << "after check" << "\n";

std::cout << "About to loop over svtxtrackseed" << std::endl;
for(auto track_iter = m_svtxtrackseed->begin(); track_iter != m_svtxtrackseed->end(); ++track_iter){
  
  const auto& trackseed = *track_iter;
  trackseed->identify();
  std::cout << std::endl;


  
  TrackSeed *tpcseedfromsvtx =  _tpc_seeds->get(trackseed->get_tpc_seed_index());
  //TrackSeed *siliconseedfromsvtx =  _tpc_seeds->get(trackseed->get_silicon_seed_index());

  if(tpcseedfromsvtx){
    std::cout << "About to identify TPC seed" << std::endl;
    //tpcseedfromsvtx->identify();

    for( auto key_iter = tpcseedfromsvtx->begin_cluster_keys(); key_iter != tpcseedfromsvtx->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      std::cout << "TPC from Svtx Cluster Key: " << cluster_key << std::endl;
      std::cout << "TrkrId: " << unsigned(TrkrDefs::getTrkrId(cluster_key)) << std::endl;    
      std::cout << "Layer: " << unsigned(TrkrDefs::getLayer(cluster_key)) << std::endl;
    }
  }
  
  /*
  if(siliconseedfromsvtx){
    std::cout << "About to identify Silicon seed" << std::endl;
    //siliconseedfromsvtx->identify();

    for( auto key_iter = siliconseedfromsvtx->begin_cluster_keys(); key_iter != siliconseedfromsvtx->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      std::cout << "Silicon from Svtx Cluster Key: " << cluster_key << std::endl;
      std::cout << "TrkrId: " << unsigned(TrkrDefs::getTrkrId(cluster_key)) << std::endl;    
      std::cout << "Layer: " << unsigned(TrkrDefs::getLayer(cluster_key)) << std::endl;    
    }
  }
  */

}
std::cout << "Finished loop over svtxtrackseed" << std::endl;
//Make TClonesArray
TClonesArray& TrkContainer = *m_container->arrTrkContainer;
TrkContainer.Clear();

//need to clear keys to make next event not
m_container->clearArrayKeys();

Int_t iTrk = 0;
  long unsigned int iKey = 0;

  //std::cout << "Before loop" << "\n";
  for( const auto& trackpair:*m_track_map )
  {
    std::cout << "start of loop" << "\n";
    //new(trkrDST[iCluster]) TrkrClusterv4();
    // TrkrClusterv4* clsArr = (TrkrClusterv4*) trkrDST.ConstructedAt(iCluster);
    

    unsigned int key = trackpair.first;
    const auto track = trackpair.second;
    //auto track_struct = create_track( track );
    /*
    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct.contributors = contributors;
    
    // get particle
    auto particle = m_g4truthinfo->GetParticle(id);
    track_struct.embed = get_embed(particle);
    add_truth_information(track_struct, particle);
    */
    

    //Store information from Track into new container that will also hold cluster information
    new(TrkContainer[iTrk]) SvtxTrackArray_v1;
    SvtxTrackArray_v1* trackContainer = (SvtxTrackArray_v1*) TrkContainer.ConstructedAt(iTrk);

    if(m_container->getArrayKeysSize()<(iKey+1)){
        m_container->resizeArrayKeys(iKey+1);
    }
    m_container->setArrayKeys(iKey,key);
    

    //Save all information from a Track into the new container
    std::cout << "track get id:" << track->get_id() <<"\n"; 
    trackContainer->set_id(track->get_id());
    trackContainer->set_vertex_id(track->get_vertex_id());
    trackContainer->set_positive_charge(track->get_positive_charge());
    trackContainer->set_chisq(track->get_chisq());
    trackContainer->set_ndf(track->get_ndf());
    trackContainer->set_crossing(track->get_crossing());

    trackContainer->set_x(track->get_x());
    trackContainer->set_y(track->get_y());
    trackContainer->set_z(track->get_z());
    trackContainer->set_px(track->get_px());
    trackContainer->set_py(track->get_py());
    trackContainer->set_pz(track->get_pz());

    TrackSeed* TPCSeed = track->get_tpc_seed();
    TrackSeed* SiliconSeed = track->get_silicon_seed();

    if(!TPCSeed){
      std::cout << "TPCSeed does not exist \n";
      trackContainer->set_does_tpc_seed_exist(false);
    }else{
    //store information from tpc seed
    trackContainer->set_does_tpc_seed_exist(true);
    std::cout << "tpcseed:" << TPCSeed->get_qOverR() <<"\n"; 
    trackContainer->tpc_seed_set_qOverR(TPCSeed->get_qOverR());
    trackContainer->tpc_seed_set_X0(TPCSeed->get_X0());
    trackContainer->tpc_seed_set_Y0(TPCSeed->get_Y0());
    trackContainer->tpc_seed_set_slope(TPCSeed->get_slope());
    trackContainer->tpc_seed_set_Z0(TPCSeed->get_Z0());
    trackContainer->tpc_seed_set_crossing(TPCSeed->get_crossing());
    }

    if(!SiliconSeed){
      std::cout << "SiliconSeed does not exist \n";
      trackContainer->set_does_silicon_seed_exist(false);
    }else{
    trackContainer->set_does_silicon_seed_exist(true);
    //store information from silicon seed
    std::cout << "silconseed:" << SiliconSeed->get_qOverR() <<"\n"; 
    trackContainer->silicon_seed_set_qOverR(SiliconSeed->get_qOverR());
    trackContainer->silicon_seed_set_X0(SiliconSeed->get_X0());
    trackContainer->silicon_seed_set_Y0(SiliconSeed->get_Y0());
    trackContainer->silicon_seed_set_slope(SiliconSeed->get_slope());
    trackContainer->silicon_seed_set_Z0(SiliconSeed->get_Z0());
    trackContainer->silicon_seed_set_crossing(SiliconSeed->get_crossing());
    }

    //trackContainer->copy_states(track);
    //store clusterkeyset from each seed
  
  //set all clusterkeys = 0
  for(int layer = 0; layer < 59; layer++){
      //trackContainer->setClusKey(layer, 0);
  }

  

  if(trackContainer->get_does_tpc_seed_exist()){
    std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
  TPCSeed->identify();
    for( auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      //store information in track array
      std::cout << "TPC clusterkey: " << cluster_key <<"\n";
      std::cout << "TPC subsurfkey: " << cluster->getSubSurfKey() << std::endl;
      uint8_t layer = TrkrDefs::getLayer(cluster_key);

      std::cout << "TPC layer: " << unsigned(TrkrDefs::getLayer(cluster_key)) << std::endl;
      //std::cout << "cluster position 0: " << cluster->getPosition(0) <<"\n";
      trackContainer->setLocalX(layer, cluster->getPosition(0));
      trackContainer->setLocalY(layer, cluster->getPosition(1));
      trackContainer->setSubSurfKey(layer, cluster->getSubSurfKey());
      trackContainer->setPhiError(layer, cluster->getRPhiError());
      trackContainer->setZError(layer, cluster->getZError());
      trackContainer->setAdc(layer, cluster->getAdc());
      trackContainer->setMaxAdc(layer, cluster->getMaxAdc());
      trackContainer->setPhiSize(layer, cluster->getPhiSize());
      trackContainer->setZSize(layer, cluster->getZSize());
      trackContainer->setOverlap(layer, cluster->getOverlap());
      trackContainer->setEdge(layer, cluster->getEdge());

      std::cout << "TPC Side: " << unsigned(TrkrDefs::getZElement(cluster_key)) << std::endl;
      std::cout << "TPC Sector Id: " << unsigned(TrkrDefs::getPhiElement(cluster_key)) << std::endl;

      trackContainer->setValid(layer, true);

      if(trackContainer->getValid(layer)){
        std::cout << "Valid is true" << std::endl;
      }

      //trackContainer->setSide(layer, TrkrDefs::getZElement(cluster_key));
      //trackContainer->setSectorId(layer, TrkrDefs::getPhiElement(cluster_key));
      //trackContainer->setClusKey(layer, cluster_key);

      //TrkrDefs::
      
    }
  }

  
  if(trackContainer->get_does_silicon_seed_exist()){
    std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
  SiliconSeed->identify();
    for( auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter )
    {
      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }
      //store cluster key
      //silicon_cluster_keys->insert_cluster_key(cluster_key);

      //store information in track array
      std::cout << "Silicon clusterkey: " << cluster_key <<"\n";

      std::cout << "Silicon subsurfkey: " << cluster->getSubSurfKey() << std::endl;
      uint8_t layer = TrkrDefs::getLayer(cluster_key);

      std::cout << "Silicon layer: " << unsigned(TrkrDefs::getLayer(cluster_key)) << std::endl;
      //std::cout << "cluster position 0: " << cluster->getPosition(0) <<"\n";
      trackContainer->setLocalX(layer, cluster->getPosition(0));
      trackContainer->setLocalY(layer, cluster->getPosition(1));
      trackContainer->setSubSurfKey(layer, cluster->getSubSurfKey());
      trackContainer->setPhiError(layer, cluster->getRPhiError());
      trackContainer->setZError(layer, cluster->getZError());
      trackContainer->setAdc(layer, cluster->getAdc());
      trackContainer->setMaxAdc(layer, cluster->getMaxAdc());
      trackContainer->setPhiSize(layer, cluster->getPhiSize());
      trackContainer->setZSize(layer, cluster->getZSize());
      trackContainer->setOverlap(layer, cluster->getOverlap());
      trackContainer->setEdge(layer, cluster->getEdge());

      std::cout << "Silicon Side: " << unsigned(TrkrDefs::getZElement(cluster_key)) << std::endl;
      std::cout << "Silicon Sector Id: " << unsigned(TrkrDefs::getPhiElement(cluster_key)) << std::endl;

      //trackContainer->setSide(layer, TrkrDefs::getZElement(cluster_key));
      //trackContainer->setSectorId(layer, TrkrDefs::getPhiElement(cluster_key));

      //trackContainer->setClusKey(layer, cluster_key);

      
    }
  }
    //Start a loop over clusters to store them in container associated with each track
    //ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const override { return m_cluster_keys.find(clusterid); }
    /*
    
    MAKE A LOOP OVER CLUSTERS HERE
    Maybe go into evaluate_clusters
    */
    
    /*
    struct cluster_with_key{
    //below is in trkclusterv5.h
    float m_local[2];          //< 2D local position [cm] 2 * 32 64bit  - cumul 1*64  
    TrkrDefs::subsurfkey m_subsurfkey; //< unique identifier for hitsetkey-surface maps 16 bit
    float m_phierr;
    float m_zerr;
    unsigned short int m_adc;           //< cluster sum adc 16
    unsigned short int m_maxadc;           //< cluster sum adc 16
    char m_phisize; // 8bit
    char m_zsize;   // 8bit
    char m_overlap; // 8bit 
    char m_edge;    // 8bit - cumul 2*64

    TrkrDefs::hitsetkey m_hitsetkey;
    };
    */
    

    

    /*
    trackContainer->clusterArray[layer].m_local[0] = cluster->getPosition(0);
    trackContainer->clusterArray[layer].m_local[1] = cluster->getPosition(1);
    trackContainer->clusterArray[layer].m_subsurfkey = cluster->getSubSurfKey();
    trackContainer->clusterArray[layer].m_phierr = cluster->getRPhiError();
    trackContainer->clusterArray[layer].m_zerr = cluster->getZError();
    trackContainer->clusterArray[layer].m_adc = cluster->getAdc();
    trackContainer->clusterArray[layer].m_maxadc = cluster->getMaxAdc();
    trackContainer->clusterArray[layer].m_phisize = cluster->getPhiSize();
    trackContainer->clusterArray[layer].m_zsize = cluster->getZSize();
    trackContainer->clusterArray[layer].m_overlap = cluster->getOverlap();
    trackContainer->clusterArray[layer].m_edge = cluster->getEdge();
    */
    


    // running iterator over track states, used to match a given cluster to a track state
    



  /*
    // loop over clusters
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      //store information in track array
      layer = TrkrDefs::getLayer(cluster_key);
      trackContainer->setLocalX(layer, cluster->getPosition(0));
      trackContainer->setLocalY(layer, cluster->getPosition(1));
      trackContainer->setSubSurfKey(layer, cluster->getSubSurfKey());
      trackContainer->setPhiError(layer, cluster->getRPhiError());
      trackContainer->setZError(layer, cluster->getZError());
      trackContainer->setAdc(layer, cluster->getAdc());
      trackContainer->setMaxAdc(layer, cluster->getMaxAdc());
      trackContainer->setPhiSize(layer, cluster->getPhiSize());
      trackContainer->setZSize(layer, cluster->getZSize());
      trackContainer->setOverlap(layer, cluster->getOverlap());
      trackContainer->setEdge(layer, cluster->getEdge());
*/
       
      /*
      // find track state that is the closest to cluster
      // this assumes that both clusters and states are sorted along r
      const auto radius( cluster_struct.r );
      float dr_min = -1;
      for( auto iter = state_iter; iter != track->end_states(); ++iter )
      {
        const auto dr = std::abs( radius - get_r( iter->second->get_x(), iter->second->get_y() ) );
        if( dr_min < 0 || dr < dr_min )
        {
          state_iter = iter;
          dr_min = dr;
        } else break;
      }

      // store track state in cluster struct
      add_trk_information( cluster_struct, state_iter->second );

      // add to track
      track_struct.clusters.push_back( cluster_struct );
      */
      ++iTrk;
      ++iKey;
      std::cout << "end of loop" << "\n";
    }
    

     //m_track_map->clear();
     m_track_map->Reset();
     m_cluster_map->Reset();
     
     //m_cluster_map->Clear();


    //cluster->getSubSurfKey();
    //get a cluskey from a subsurfkey
    //TPCSeed->find_cluster_key(cluskey);

    



}

//_____________________________________________________________________
void DSTWriter::evaluate_event()
{
  if(!m_container) return;

  // create event struct
  DSTContainerv3::EventStruct event;
  if( m_hitsetcontainer )
  {
    // loop over hitsets
    for(const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
    {
      const auto trkrId = TrkrDefs::getTrkrId(hitsetkey);
      const auto layer = TrkrDefs::getLayer(hitsetkey);
      assert(layer<DSTContainerv3::EventStruct::max_layer);

      if(m_cluster_map)
      {

        // fill cluster related information
        const auto clusters = m_cluster_map->getClusters(hitsetkey);
        const int nclusters = std::distance( clusters.first, clusters.second );
        
        switch( trkrId )
        {
          case TrkrDefs::mvtxId: event.nclusters_mvtx += nclusters; break;
          case TrkrDefs::inttId: event.nclusters_intt += nclusters; break;
          case TrkrDefs::tpcId: event.nclusters_tpc += nclusters; break;
          case TrkrDefs::micromegasId: event.nclusters_micromegas += nclusters; break;
        }

        event.nclusters[layer] += nclusters;
      }
    }
  }

  // store
  m_container->addEvent(event);
}

//_____________________________________________________________________
void DSTWriter::evaluate_clusters()
{
  std::cout << "entering clusters" << std::endl;

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
  


  // clear array
  m_container->clearClusters();
  //see if we need to reset this new map
  m_cluster_map_arr->Reset();

  // TClonesArray& ar = *arrCls;
  // ar.Clear();

  //TClonesArray& arDST = *m_container->arrClsDST;
  //arDST.Clear();

  //TClonesArray& trkrDST = *m_container->trkrClsDST;
  //trkrDST.Clear();

  //TClonesArray& arrKeyDST = *m_container->arrKeyDST;
  //arrKeyDST.Clear();

  std::cout << "DSTWriter::evaluate_clusters  current clusters: " << m_cluster_map->size() << std::endl;
  Int_t iCluster = 0;
  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    std::cout << "looping over hitsetkey " << std::hex << (unsigned) hitsetkey << "\n";
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // create cluster structure
      // auto cluster_struct = create_cluster( key, cluster );
      // add_cluster_size( cluster_struct, key, m_cluster_hit_map );
      // add_cluster_energy( cluster_struct, key, m_cluster_hit_map, m_hitsetcontainer );
      // truth information
      // const auto g4hits = find_g4hits( key );
      // add_truth_information( cluster_struct, g4hits );

      // add in array
      // m_container->addCluster( cluster_struct );
      // std::cout << "added cluster with key " << std::hex << (ulong) key << "\n";

      // new(ar[iCluster]) DSTContainerTcl::ClusterStruct( key, cluster );
      // new(arDST[iCluster]) DSTContainerv3::ClusterStruct( key, cluster );
      //TrkrClusterv4 clusIn;
      //clusIn.CopyFrom(cluster);

      // TrkrClusterv4 *clusIn = dynamic_cast<TrkrClusterv4*> (cluster);
      // new(arDST[iCluster]) DSTContainerv3::ClusterStruct( hitsetkey, clusIn );
      // new(trkrDST[iCluster]) TrkrClusterv4();

      // hitsetkey
      // new(arrKeyDST[iCluster]) TrkrDefs::hitsetkey(hitsetkey);
      // new(arrKeyDST[iCluster]) DSTContainerv3::ClusterKeyStruct(hitsetkey);
      // new(arDST[iCluster]) DSTContainerv3::ClusterStruct(hitsetkey, clusIn);
      //new(arDST[iCluster]) DSTContainerv3::ClusterStruct(hitsetkey, clusIn, key);
      // TrkrClusterv4* clsArr = (TrkrClusterv4*) trkrDST.ConstructedAt(iCluster);
      // clsArr->CopyFrom(cluster);
      std::cout << "adding cluster to new container" << std::endl;
      //std::cout << key << std::endl;
      m_cluster_map_arr->addClusterSpecifyKey(key, cluster);
      ++iCluster;
    }
  }
  // tcl->Fill();
  //reset cluster map
  m_cluster_map->Reset();
  //std::cout << "DSTWriter::evaluate_clusters - clusters: " << m_container->clusters().size() << std::endl;
  
}

//_____________________________________________________________________
/* Don't need this because it will all go in m_container
void DSTWriter::evaluate_trackarray()
{
  //new setup for trackarray
  if(!(m_track_array && m_container)) return;

  //clear array
  m_container->clearTracks();
  //
  TClonesArray& TrkArr = *m_container->arrTrkArr;
  TrkArr.Clear();

}
*/
//_____________________________________________________________________
void DSTWriter::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  //clear TrackArray
  TClonesArray& TrkArr = *m_container->arrTrkArr;
  TrkArr.Clear();

  //Clear our new Track container
  //TClonesArray& TrkContainer = *m_container->TrkArrayContainer;
  //TrkContainer.Clear();

  //std::vector<unsigned int>& ArrayKeys = *m_container->TrackArrayKeys;
  m_container->clearArrayKeys();
  // clear array
  m_container->clearTracks();

  // TClonesArray& ar = *arrTrk;
  // ar.Clear();

  // TClonesArray& arDST = *m_container->arrClsDST;
  // arDST.Clear();

  Int_t iTrk = 0;
  long unsigned int iKey = 0;
  for( const auto& trackpair:*m_track_map )
  {

    //new(trkrDST[iCluster]) TrkrClusterv4();
    // TrkrClusterv4* clsArr = (TrkrClusterv4*) trkrDST.ConstructedAt(iCluster);
    

    unsigned int key = trackpair.first;
    const auto track = trackpair.second;
    auto track_struct = create_track( track );
    /*
    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct.contributors = contributors;
    
    // get particle
    auto particle = m_g4truthinfo->GetParticle(id);
    track_struct.embed = get_embed(particle);
    add_truth_information(track_struct, particle);
    */
    /*
    // running iterator over track states, used to match a given cluster to a track state
    auto state_iter = track->begin_states();

    // loop over clusters
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "DSTWriter::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      // create new cluster struct
      auto cluster_struct = create_cluster( cluster_key, cluster );
      add_cluster_size( cluster_struct, cluster_key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, cluster_key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( cluster_key );
      add_truth_information( cluster_struct, g4hits );

      // find track state that is the closest to cluster
      // this assumes that both clusters and states are sorted along r
      const auto radius( cluster_struct.r );
      float dr_min = -1;
      for( auto iter = state_iter; iter != track->end_states(); ++iter )
      {
        const auto dr = std::abs( radius - get_r( iter->second->get_x(), iter->second->get_y() ) );
        if( dr_min < 0 || dr < dr_min )
        {
          state_iter = iter;
          dr_min = dr;
        } else break;
      }

      // store track state in cluster struct
      add_trk_information( cluster_struct, state_iter->second );

      // add to track
      track_struct.clusters.push_back( cluster_struct );
    }
    */

    //Store information from Track into new container that will also hold cluster information
    /*
    new(TrkContainer[iTrk]) SvtxTrackArray_v1;
    SvtxTrackArray_v1* trackContainer = (SvtxTrackArray_v1*) TrkContainer.ConstructedAt(iTrk);


    //Save all information from a Track into the new container
     
    trackContainer->set_id(track->get_id());
    trackContainer->set_vertex_id(track->get_vertex_id);
    trackContainer->set_positive_charge(track->get_positive_charge);
    trackContainer->set_chisq(track->get_chisq());
    trackContainer->set_ndf(track->get_ndf());
    trackContainer->set_crossing(track->get_crossing());

    TrackSeed* TPCSeed = track.get_tpc_seed;
    TrackSeed* SiliconSeed = track->get_silicon_seed;
    */ 

    /*
    ClusterKeySet m_cluster_keys;
  
    float m_qOverR = NAN;
    float m_X0 = NAN;
    float m_Y0 = NAN;
    float m_slope = NAN;
    float m_Z0 = NAN;

    short int m_crossing = std::numeric_limits<short int>::max();
    */

    /*
    //store information from tpc seed
    trackContainer->tpc_seed_set_qOverR(TPCSeed->get_qOverR());
    trackContainer->tpc_seed_set_X0(TPCSeed->get_X0());
    trackContainer->tpc_seed_set_Y0(TPCSeed->get_Y0());
    trackContainer->tpc_seed_set_slope(TPCSeed->get_slope());
    trackContainer->tpc_seed_set_Z0(TPCSeed->get_Z0());
    trackContainer->tpc_seed_set_crossing(TPCSeed->get_crossing());

    //store information from silicon seed
    trackContainer->silicon_seed_set_qOverR(SiliconSeed->get_qOverR());
    trackContainer->silicon_seed_set_X0(SiliconSeed->get_X0());
    trackContainer->silicon_seed_set_Y0(SiliconSeed->get_Y0());
    trackContainer->silicon_seed_set_slope(SiliconSeed->get_slope());
    trackContainer->silicon_seed_set_Z0(SiliconSeed->get_Z0());
    trackContainer->silicon_seed_set_crossing(SiliconSeed->get_crossing());
    */

    //Start a loop over clusters to store them in container associated with each track
    //ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const override { return m_cluster_keys.find(clusterid); }
    /*
    
    MAKE A LOOP OVER CLUSTERS HERE
    Maybe go into evaluate_clusters
    */

    //cluster->getSubSurfKey();
    //get a cluskey from a subsurfkey
    //TPCSeed->find_cluster_key(cluskey);



    
    //m_container->addTrack( track_struct );
    // new(ar[iTrk]) DSTContainerTcl::TrackStruct(track);

    new(TrkArr[iTrk]) SvtxTrack_v4;
    //SvtxTrackArray_v1();
    SvtxTrack_v4* trackArr = (SvtxTrack_v4*) TrkArr.ConstructedAt(iTrk);
    trackArr->CopyFrom(trackpair.second);
 
    // if(ArrayKeys.size()<iKey){
    //    ArrayKeys.resize(iKey);
   // }
   // ArrayKeys[iKey] = key;
    if(m_container->getArrayKeysSize()<(iKey+1)){
        m_container->resizeArrayKeys(iKey+1);
    }
    m_container->setArrayKeys(iKey,key);


    ++iTrk;
    ++iKey;
  }
  //Clear Track Map to free up space
  m_track_map->clear();
  
  std::cout << "DSTWriter::evaluate_tracks - tracks: " << m_container->tracks().size() << std::endl;

  
}

//_____________________________________________________________________
DSTWriter::G4HitSet DSTWriter::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  // check maps
  if( !( m_cluster_hit_map && m_hit_truth_map ) ) return G4HitSet();

  // check if in map
  auto map_iter = m_g4hit_map.lower_bound( cluster_key );
  if( map_iter != m_g4hit_map.end() && cluster_key == map_iter->first )
  { return map_iter->second; }

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
  for(const auto& [first,hit_key]:range_adaptor( m_cluster_hit_map->getHits(cluster_key)))
  {

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit* g4hit = nullptr;

      switch( TrkrDefs::getTrkrId( hitset_key ) )
      {
        case TrkrDefs::mvtxId:
        if( m_g4hits_mvtx ) g4hit = m_g4hits_mvtx->findHit( g4hit_key );
        break;

        case TrkrDefs::inttId:
        if( m_g4hits_intt ) g4hit = m_g4hits_intt->findHit( g4hit_key );
        break;

        case TrkrDefs::tpcId:
        if( m_g4hits_tpc ) g4hit = m_g4hits_tpc->findHit( g4hit_key );
        break;

        case TrkrDefs::micromegasId:
        if( m_g4hits_micromegas ) g4hit = m_g4hits_micromegas->findHit( g4hit_key );
        break;

        default: break;
      }

      if( g4hit ) out.insert( g4hit );
      else std::cout << "DSTWriter::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> DSTWriter::get_max_contributor( SvtxTrack* track ) const
{
  if(!(m_track_map && m_cluster_map && m_g4truthinfo)) return {0,0};

  // maps MC track id and number of matching g4hits
  using IdMap = std::map<int,int>;
  IdMap contributor_map;

  // loop over clusters
  for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
  {
    const auto& cluster_key = *key_iter;
    for( const auto& hit:find_g4hits( cluster_key ) )
    {
      const int trkid = hit->get_trkid();
      auto iter = contributor_map.lower_bound( trkid );
      if( iter == contributor_map.end() || iter->first != trkid )
      {
        contributor_map.insert(iter, std::make_pair(trkid,1));
      } else ++iter->second;
    }
  }

  if( contributor_map.empty() ) return {0,0};
  else return *std::max_element(
    contributor_map.cbegin(), contributor_map.cend(),
    []( const IdMap::value_type& first, const IdMap::value_type& second )
    { return first.second < second.second; } );

}

//_____________________________________________________________________
int DSTWriter::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

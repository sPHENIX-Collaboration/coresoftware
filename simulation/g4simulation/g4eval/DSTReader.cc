/*!
 * \file DSTReader.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "DSTReader.h"
#include "DSTContainerv1.h"

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
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>

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
DSTReader::DSTReader( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int DSTReader::Init(PHCompositeNode* topNode )
{

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

  // auto newNode = new PHIODataNode<PHObject>( new DSTContainerv1, "DSTContainer","PHObject");
  // evalNode->addNode(newNode);

  // // DST container
  // m_container = findNode::getClass<DSTContainerv1>(topNode, "DSTContainer");


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
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output container
  // if( m_container ) m_container->Reset();

  // if(m_flags&WriteEvent) 
  //   evaluate_event();
  // if(m_flags&WriteClusters) 
    evaluate_clusters();
  // if(m_flags&WriteTracks)
    // evaluate_tracks();
  std::cout << "evaluate cluster finished" << "\n";


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
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap_2");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<DSTContainerv1>(topNode, "DSTContainer");

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void DSTReader::evaluate_event()
{
  if(!m_container) return;

  // create event struct

  // store

  // m_container->addEvent(event);
}

//_____________________________________________________________________
void DSTReader::evaluate_clusters()
{

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  // m_container->clearClusters();

  std::cout << "DSTReader::evaluate_clusters - clusters: " << m_container->clusters().size() << std::endl;
  // debugging
  showAll();
  m_cluster_map->Reset();
  std::cout << "DSTReader::evaluate_clusters - cleared. current clusters: " << m_cluster_map->size() << std::endl;

  uint8_t id = 0;
  for (auto& clusterStruct : m_container->clusters())
  {
    std::cout << "cluster " << (unsigned) id << "\n";
    // TrkrCluster newCluster = recover_cluster(cluster);

    auto newCluster = new TrkrClusterv3;
    newCluster->setLocalX(clusterStruct.loc_x);
    newCluster->setLocalY(clusterStruct.loc_y);
    TrkrDefs::hitsetkey key = 0;
    // std::cout << std::bitset<32>(key) << "\n";
    // unsigned char z_unsigned = (unsigned) clusterStruct.z_seg;
    // key |= (z_unsigned << TrkrDefs::kBitShiftZElement);
    key |= (clusterStruct.z_seg << TrkrDefs::kBitShiftZElement);
    // std::cout << std::bitset<32>(key) << "\n";
    // unsigned char phi_unsigned = (unsigned) clusterStruct.phi_seg;
    // key |= (phi_unsigned << TrkrDefs::kBitShiftPhiElement);
    key |= (clusterStruct.phi_seg << TrkrDefs::kBitShiftPhiElement);
    // std::cout << std::bitset<32>(key) << "\n";
    // unsigned char layer_unsigned = (unsigned) clusterStruct.layer;
    // key |= (layer_unsigned << TrkrDefs::kBitShiftLayer);
    key |= (clusterStruct.layer << TrkrDefs::kBitShiftLayer);
    // std::cout << std::bitset<32>(key) << "\n";
    key |= (id << TrkrDefs::kBitShiftTrkrId);
    // std::cout << std::bitset<32>(key) << "\n";
    TrkrDefs::cluskey cluskey = key;
    cluskey = (cluskey << TrkrDefs::kBitShiftClusId);
    // std::cout << std::bitset<64>(cluskey) << "\n";


    // set size and error
    // for (int j = 0; j < 3; ++j) {
    //   for (int i = 0; i < 3; ++i) {
    //     newCluster->setSize(i, j, DSTContainerv1::covarIndex(i, j));
    //     newCluster->setError(i, j, DSTContainerv1::covarIndex(i, j));
    //   }
    // }

    int nLocal = 2;
    for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
        newCluster->setActsLocalError(iLocal, jLocal,
                                     clusterStruct.actsLocalError[iLocal][jLocal]);
        // std::cout << "actslocalerror:" << newCluster->
          // getActsLocalError(iLocal, jLocal) << "\n";
      }
    }
    newCluster->setSubSurfKey(clusterStruct.subSurfKey);
    // std::cout << "subsurfkey: " << newCluster->getSubSurfKey() << "\n";

    newCluster->setAdc(clusterStruct.adc);
    // std::cout << "adc: " << newCluster->getAdc() << "\n";


    // debugging
    // std::cout << "key from clust: " << std::hex << clusterStruct.clusterKey << "\n";
    // std::cout << "key by me: " << std::hex << cluskey << std::dec << "\n";
    // for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
    //   std::cout << "positions: " << newCluster->getPosition(iLocal) << "\n";
    // }
    // std::cout << "rphierror: " << newCluster->getRPhiError() << "\n";
    // std::cout << "zerror: " << newCluster->getZError() << "\n";





    // get hitsetkey from cluster
    // const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( cluskey );
    // // get cluster index in vector
    // const auto index = TrkrDefs::getClusIndex( cluskey );
    // std::cout << "index:" << index << "\n";

    // m_cluster_map->addClusterSpecifyKey(cluskey, &newCluster);
    m_cluster_map->addClusterSpecifyKey(clusterStruct.clusterKey, newCluster);
    ++id;
    std::cout << "DSTReader::evaluate_clusters - saved clusters: " << m_cluster_map->size() << std::endl;
    showMe();
  }

  std::cout << "DSTReader::evaluate_clusters - saved clusters: " << m_cluster_map->size() << std::endl;

  showAll();

}

//_____________________________________________________________________
void DSTReader::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  std::cout << "reading file" << "\n";

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
  std::cout << "DSTReader::evaluate_tracks - saved tracks: " << m_track_map->size() << std::endl;


}

//_____________________________________________________________________
DSTReader::G4HitSet DSTReader::find_g4hits( TrkrDefs::cluskey cluster_key ) const
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
      else std::cout << "DSTReader::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> DSTReader::get_max_contributor( SvtxTrack* track ) const
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

// //! create svx track from struct
// SvtxTrack DSTReader::recover_track(DSTContainerv1::TrackStruct trackStruct)
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

// TrkrCluster DSTReader::recover_cluster(DSTContainerv1::ClusterStruct clusterStruct)
// {
//   TrkrClusterv3 cluster;
//   cluster.setLocalX(clusterStruct.loc_x);
//   cluster.setLocalY(clusterStruct.loc_y);

//   return cluster;
// }

//_____________________________________________________________________
int DSTReader::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }


void DSTReader::printCluster(TrkrCluster& newCluster) const {
    int nLocal = 2;
    for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
        std::cout << "actslocalerror:" << newCluster.
          getActsLocalError(iLocal, jLocal) << std::endl;
      }
    }
    std::cout << "subsurfkey: " << newCluster.getSubSurfKey() << std::endl;
    std::cout << "adc: " << newCluster.getAdc() << std::endl;
    for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
      std::cout << "positions: " << newCluster.getPosition(iLocal) << std::endl;
    }
    std::cout << "rphierror: " << newCluster.getRPhiError() << std::endl;
    std::cout << "zerror: " << newCluster.getZError() << std::endl;
}

void DSTReader::showMe() const {
  std::cout << "hitset size: " << m_hitsetcontainer->size() << std::endl;
  for( const auto& [hitsetkey, hitset]:range_adaptor(m_hitsetcontainer->getHitSets())) {
    std::cout << "looping over hitsetkey " << std::hex << hitsetkey << std::endl;
    for( const auto& [key, cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey))) {
      std::cout << "printCluster " << std::hex << key << std::dec << std::endl;
      // printCluster(*cluster);
      // std::cout << "cluster is valid " << cluster->isValid() << std::endl;
      cluster->identify();
    }
    break;
  }
}


void DSTReader::showAll() const {
  std::cout << "show all clusters of " << m_cluster_map->size() << std::endl;
  for (const auto& hitset : m_cluster_map->getHitSetKeys()) {
    for( const auto& [key, cluster]:range_adaptor(m_cluster_map->getClusters(hitset))) {
      std::cout << "identify cluster with key " << std::hex << key << std::dec << std::endl;
      std::cout << "cluster is valid? " << cluster->isValid() << std::endl;
      cluster->identify();
      std::cout << "our time" << "\n";
      printCluster(*cluster);
    }
  }
}

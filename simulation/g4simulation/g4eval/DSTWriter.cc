/*!
 * \file DSTWriter.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "DSTWriter.h"
#include "DSTContainerv3.h"
#include "DSTContainerTcl.h"

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
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

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

  // TClonesArary
  // fcl = new TFile("dstcl.root", "recreate");
  // tcl = new TTree("tcl", "dst tree");
  // arrEvt = new TClonesArray("DSTContainerv3::EventStruct");
  // arrTrk = new TClonesArray("DSTContainerv3::TrackStruct");
  // arrCls = new TClonesArray("DSTContainerv3::ClusterStruct");
  // tcl->Branch("evt", &arrEvt);
  // tcl->Branch("trk", &arrTrk);
  // tcl->Branch("cls", &arrCls);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTWriter::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTWriter::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output container
  if( m_container ) m_container->Reset();
  std::cout << "DSTWriter::process_event" << std::endl;

  if(m_flags&WriteEvent)
    evaluate_event();
  if(m_flags&WriteClusters) {
    std::cout << "DSTWriter::process_event doing clusters" << std::endl;
    evaluate_clusters();
  }
  if(m_flags&WriteTracks) {
    std::cout << "DSTWriter::process_event doing tracks" << std::endl;
    // evaluate_tracks();
  }
  std::cout << "exiting event" << "\n";


  // clear maps
  m_g4hit_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTWriter::End(PHCompositeNode* )
{

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

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<DSTContainerv3>(topNode, "DSTContainer");

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

  // TClonesArray& ar = *arrCls;
  // ar.Clear();

  TClonesArray& arDST = *m_container->arrClsDST;
  arDST.Clear();

  TClonesArray& trkrDST = *m_container->trkrClsDST;
  trkrDST.Clear();

  TClonesArray& arrKeyDST = *m_container->arrKeyDST;
  arrKeyDST.Clear();

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
      TrkrClusterv4 clusIn;
      clusIn.CopyFrom(cluster);

      // TrkrClusterv4 *clusIn = dynamic_cast<TrkrClusterv4*> (cluster);
      // new(arDST[iCluster]) DSTContainerv3::ClusterStruct( hitsetkey, clusIn );
      // new(trkrDST[iCluster]) TrkrClusterv4();

      // hitsetkey
      // new(arrKeyDST[iCluster]) TrkrDefs::hitsetkey(hitsetkey);
      // new(arrKeyDST[iCluster]) DSTContainerv3::ClusterKeyStruct(hitsetkey);
      // new(arDST[iCluster]) DSTContainerv3::ClusterStruct(hitsetkey, clusIn);
      new(arDST[iCluster]) DSTContainerv3::ClusterStruct(hitsetkey, clusIn, key);
      // TrkrClusterv4* clsArr = (TrkrClusterv4*) trkrDST.ConstructedAt(iCluster);
      // clsArr->CopyFrom(cluster);
      ++iCluster;
    }
  }
  // tcl->Fill();

  std::cout << "DSTWriter::evaluate_clusters - clusters: " << m_container->clusters().size() << std::endl;
  
}

//_____________________________________________________________________
void DSTWriter::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;


  // clear array
  m_container->clearTracks();

  // TClonesArray& ar = *arrTrk;
  // ar.Clear();

  // TClonesArray& arDST = *m_container->arrClsDST;
  // arDST.Clear();

  Int_t iTrk = 0;
  for( const auto& trackpair:*m_track_map )
  {


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
    m_container->addTrack( track_struct );
    // new(ar[iTrk]) DSTContainerTcl::TrackStruct(track);
    ++iTrk;
  }
  
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

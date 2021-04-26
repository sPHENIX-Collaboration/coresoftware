/*!
 * \file TrackEvaluation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluation.h"
#include "TrackEvaluationContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <intt/InttDefs.h>
#include <micromegas/MicromegasDefs.h>
#include <mvtx/MvtxDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <tpc/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
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
  TrackEvaluationContainerv1::TrackStruct create_track( SvtxTrack* track )
  {
    TrackEvaluationContainerv1::TrackStruct trackStruct;

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
    trackStruct.r = get_r( trackStruct.x, trackStruct.y );
    trackStruct.phi = std::atan2( trackStruct.y, trackStruct.x );

    trackStruct.px = track->get_px();
    trackStruct.py = track->get_py();
    trackStruct.pz = track->get_pz();
    trackStruct.pt = get_pt( trackStruct.px, trackStruct.py );
    trackStruct.p = get_p( trackStruct.px, trackStruct.py, trackStruct.pz );
    trackStruct.eta = get_eta( trackStruct.p, trackStruct.pz );

    return trackStruct;
  }

  //! create cluster struct from svx cluster
  TrackEvaluationContainerv1::ClusterStruct create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster )
  {
    TrackEvaluationContainerv1::ClusterStruct cluster_struct;
    cluster_struct.layer = TrkrDefs::getLayer(key);
    cluster_struct.x = cluster->getX();
    cluster_struct.y = cluster->getY();
    cluster_struct.z = cluster->getZ();
    cluster_struct.r = get_r( cluster_struct.x, cluster_struct.y );
    cluster_struct.phi = std::atan2( cluster_struct.y, cluster_struct.x );
    cluster_struct.phi_error = cluster->getPhiError();
    cluster_struct.z_error = cluster->getZError();

    return cluster_struct;
  }

  //! add track information
  void add_trk_information( TrackEvaluationContainerv1::ClusterStruct& cluster, SvtxTrackState* state )
  {
    // need to extrapolate the track state to the right cluster radius to get proper residuals  
    const auto trk_r = get_r( state->get_x(), state->get_y() );
    const auto dr = cluster.r - trk_r;
    const auto trk_drdt = (state->get_x()*state->get_px() + state->get_y()*state->get_py())/trk_r;
    const auto trk_dxdr = state->get_px()/trk_drdt;
    const auto trk_dydr = state->get_py()/trk_drdt;
    const auto trk_dzdr = state->get_pz()/trk_drdt;

    // store state position
    cluster.trk_x = state->get_x() + dr*trk_dxdr;
    cluster.trk_y = state->get_y() + dr*trk_dydr;
    cluster.trk_z = state->get_z() + dr*trk_dzdr;
    cluster.trk_r = get_r( cluster.trk_x, cluster.trk_y );
    cluster.trk_phi = std::atan2( cluster.trk_y, cluster.trk_x );

    /* store local momentum information */
    cluster.trk_px = state->get_px();
    cluster.trk_py = state->get_py();
    cluster.trk_pz = state->get_pz();

    /*
    store state angles in (r,phi) and (r,z) plans
    they are needed to study space charge distortions
    */
    const auto cosphi( std::cos( cluster.trk_phi ) );
    const auto sinphi( std::sin( cluster.trk_phi ) );
    const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto trk_pz = state->get_pz();
    cluster.trk_alpha = std::atan2( trk_pphi, trk_pr );
    cluster.trk_beta = std::atan2( trk_pz, trk_pr );
    cluster.trk_phi_error = state->get_phi_error();
    cluster.trk_z_error = state->get_z_error();

  }

  //! number of hits associated to cluster
  void add_cluster_size( TrackEvaluationContainerv1::ClusterStruct& cluster, TrkrDefs::cluskey clus_key, TrkrClusterHitAssoc* cluster_hit_map )
  {
    if( !cluster_hit_map ) return;
    const auto range = cluster_hit_map->getHits(clus_key);

    // store full size
    cluster.size =  std::distance( range.first, range.second );

    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId == TrkrDefs::micromegasId)
    {

      // for micromegas the directional cluster size depends on segmentation type
      auto segmentation_type = MicromegasDefs::getSegmentationType(clus_key);
      if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_Z ) cluster.z_size = cluster.size;
      else cluster.phi_size = cluster.size;

    } else {

      // for other detectors, one must loop over the constituting hits
      std::set<int> phibins;
      std::set<int> zbins;
      for(const auto& [first, hit_key]:range_adaptor(range))
      {
        switch( detId )
        {
          default: break;
          case TrkrDefs::mvtxId:
          {
            phibins.insert( MvtxDefs::getRow( hit_key ) );
            zbins.insert( MvtxDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::inttId:
          {
            phibins.insert( InttDefs::getRow( hit_key ) );
            zbins.insert( InttDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::tpcId:
          {
            phibins.insert( TpcDefs::getPad( hit_key ) );
            zbins.insert( TpcDefs::getTBin( hit_key ) );
            break;
          }
        }
      }
      cluster.phi_size = phibins.size();
      cluster.z_size = zbins.size();
    }
  }


  //! hit energy for a given cluster
  void add_cluster_energy( TrackEvaluationContainerv1::ClusterStruct& cluster, TrkrDefs::cluskey clus_key,
    TrkrClusterHitAssoc* cluster_hit_map,
    TrkrHitSetContainer* hitsetcontainer )
  {

    // check container
    if(!(cluster_hit_map && hitsetcontainer)) return;

    // for now this is only filled for micromegas
    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId != TrkrDefs::micromegasId) return;

    const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(clus_key);
    const auto hitset = hitsetcontainer->findHitSet( hitset_key );
    if( !hitset ) return;

    const auto range = cluster_hit_map->getHits(clus_key);
    cluster.energy_max = 0;
    cluster.energy_sum = 0;

    for( const auto& pair:range_adaptor(range))
    {
      const auto hit = hitset->getHit( pair.second );
      if( hit )
      {
        const auto energy = hit->getEnergy();
        cluster.energy_sum += energy;
        if( energy > cluster.energy_max ) cluster.energy_max = energy;
      }
    }

  }

  // add truth information
  void add_truth_information( TrackEvaluationContainerv1::ClusterStruct& cluster, std::set<PHG4Hit*> hits )
  {
    // need to extrapolate g4hits position to the cluster r
    /* this is done by using a linear extrapolation, with a straight line going through all provided G4Hits in and out positions */
    const auto rextrap = cluster.r;
    cluster.truth_size = hits.size();
    cluster.truth_x = interpolate<&PHG4Hit::get_x>( hits, rextrap );
    cluster.truth_y = interpolate<&PHG4Hit::get_y>( hits, rextrap );
    cluster.truth_z = interpolate<&PHG4Hit::get_z>( hits, rextrap );
    cluster.truth_r = get_r( cluster.truth_x, cluster.truth_y );
    cluster.truth_phi = std::atan2( cluster.truth_y, cluster.truth_x );

    /* add truth momentum information */
    cluster.truth_px = interpolate<&PHG4Hit::get_px>( hits, rextrap );
    cluster.truth_py = interpolate<&PHG4Hit::get_py>( hits, rextrap );
    cluster.truth_pz = interpolate<&PHG4Hit::get_pz>( hits, rextrap );

    /*
    store state angles in (r,phi) and (r,z) plans
    they are needed to study space charge distortions
    */
    const auto cosphi( std::cos( cluster.truth_phi ) );
    const auto sinphi( std::sin( cluster.truth_phi ) );
    const auto truth_pphi = -cluster.truth_px*sinphi + cluster.truth_py*cosphi;
    const auto truth_pr = cluster.truth_px*cosphi + cluster.truth_py*sinphi;

    cluster.truth_alpha = std::atan2( truth_pphi, truth_pr );
    cluster.truth_beta = std::atan2( cluster.truth_pz, truth_pr );
    if(std::isnan(cluster.truth_alpha) || std::isnan(cluster.truth_beta))
    {
      // recalculate
      double truth_alpha = 0;
      double truth_beta = 0;
      double sum_w = 0;
      for( const auto& hit:hits )
      {
        const auto px = hit->get_x(1) - hit->get_x(0);
        const auto py = hit->get_y(1) - hit->get_y(0);
        const auto pz = hit->get_z(1) - hit->get_z(0);
        const auto pphi = -px*sinphi + py*cosphi;
        const auto pr = px*cosphi + py*sinphi;

        const auto w =  hit->get_edep();
        if( w < 0 ) continue;

        sum_w += w;
        truth_alpha += w*std::atan2( pphi, pr );
        truth_beta += w*std::atan2( pz, pr );
      }
      truth_alpha /= sum_w;
      truth_beta /= sum_w;
      cluster.truth_alpha = truth_alpha;
      cluster.truth_beta = truth_beta;
    }

  }

  // ad}d truth information
  void add_truth_information( TrackEvaluationContainerv1::TrackStruct& track, PHG4Particle* particle )
  {
    if( particle )
    {
      track.is_primary = is_primary( particle );
      track.pid = particle->get_pid();
      track.truth_px = particle->get_px();
      track.truth_py = particle->get_py();
      track.truth_pz = particle->get_pz();
      track.truth_pt = get_pt( track.truth_px, track.truth_py );
      track.truth_p = get_p( track.truth_px, track.truth_py, track.truth_pz );
      track.truth_eta = get_eta( track.truth_p, track.truth_pz );
    }
  }

}

//_____________________________________________________________________
TrackEvaluation::TrackEvaluation( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int TrackEvaluation::Init(PHCompositeNode* topNode )
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TrackEvaluation::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "TrackEvaluation::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto newNode = new PHIODataNode<PHObject>( new TrackEvaluationContainerv1, "TrackEvaluationContainer","PHObject");
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackEvaluation::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackEvaluation::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output container
  if( m_container ) m_container->Reset();

  if(m_flags&EvalEvent) evaluate_event();
  if(m_flags&EvalClusters) evaluate_clusters();
  if(m_flags&EvalTracks) evaluate_tracks();

  // clear maps
  m_g4hit_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackEvaluation::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackEvaluation::load_nodes( PHCompositeNode* topNode )
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
  m_container = findNode::getClass<TrackEvaluationContainerv1>(topNode, "TrackEvaluationContainer");

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
void TrackEvaluation::evaluate_event()
{
  if(!m_container) return;

  // create event struct
  TrackEvaluationContainerv1::EventStruct event;
  if( m_hitsetcontainer )
  {
    // loop over hitsets
    for(const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
    {
      const auto trkrId = TrkrDefs::getTrkrId(hitsetkey);
      const auto layer = TrkrDefs::getLayer(hitsetkey);
      assert(layer<TrackEvaluationContainerv1::EventStruct::max_layer);

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
void TrackEvaluation::evaluate_clusters()
{

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  m_container->clearClusters();

  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // create cluster structure
      auto cluster_struct = create_cluster( key, cluster );
      add_cluster_size( cluster_struct, key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( key );
      add_truth_information( cluster_struct, g4hits );

      // add in array
      m_container->addCluster( cluster_struct );
    }
  }
  
  std::cout << "TrackEvaluation::evaluate_clusters - clusters: " << m_container->clusters().size() << std::endl;
  
}

//_____________________________________________________________________
void TrackEvaluation::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  // clear array
  m_container->clearTracks();

  for( const auto& trackpair:*m_track_map )
  {

    const auto track = trackpair.second;
    auto track_struct = create_track( track );

    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct.contributors = contributors;
    
    // get particle
    auto particle = m_g4truthinfo->GetParticle(id);
    track_struct.embed = get_embed(particle);
    add_truth_information(track_struct, particle);

    // running iterator over track states, used to match a given cluster to a track state
    auto state_iter = track->begin_states();

    // loop over clusters
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackEvaluation::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
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
      /* this assumes that both clusters and states are sorted along r */
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
    m_container->addTrack( track_struct );
  }
  
  std::cout << "TrackEvaluation::evaluate_tracks - tracks: " << m_container->tracks().size() << std::endl;

  
}

//_____________________________________________________________________
TrackEvaluation::G4HitSet TrackEvaluation::find_g4hits( TrkrDefs::cluskey cluster_key ) const
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
      else std::cout << "TrackEvaluation::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> TrackEvaluation::get_max_contributor( SvtxTrack* track ) const
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
int TrackEvaluation::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

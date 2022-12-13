/*!
 * \file TrackEvaluation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluation.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TVector3.h>

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
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  //! radius
  template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  //! needed for weighted linear interpolation
  struct interpolation_data_t
  {
    using list = std::vector<interpolation_data_t>;
    double x() const { return position.x(); }
    double y() const { return position.y(); }
    double z() const { return position.z(); }

    double px() const { return momentum.x(); }
    double py() const { return momentum.y(); }
    double pz() const { return momentum.z(); }

    TVector3 position;
    TVector3 momentum;
    double weight = 1;
  };

  //! calculate the average of member function called on all members in collection
  template<double (interpolation_data_t::*accessor)() const>
  double average( const interpolation_data_t::list& hits )
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      if(std::isnan(x)) continue;

      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      sw += w;
      swx += w*x;
    }

    if( !valid ) return NAN;
    return swx/sw;
  }

  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys( SvtxTrack* track )
  {
    std::vector<TrkrDefs::cluskey> out;
    for( const auto& seed: { track->get_silicon_seed(), track->get_tpc_seed() } )
    {
      if( seed )
      { std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter( out ) ); }
    }
    
    return out;
  }

  //! true if a track is a primary
  inline int is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }

  //! get mask from track clusters
  int64_t get_mask( SvtxTrack* track )
  { 
    const auto cluster_keys = get_cluster_keys( track );
    return std::accumulate( cluster_keys.begin(), cluster_keys.end(), int64_t(0),
      []( int64_t value, const TrkrDefs::cluskey& key ) {
        return TrkrDefs::getLayer(key)<64 ? value|(1LL<<TrkrDefs::getLayer(key)) : value;
      } );
  }

  //! return number of clusters of a given type
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    const auto cluster_keys = get_cluster_keys( track );
    return std::count_if( cluster_keys.begin(), cluster_keys.end(),
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

  //! number of hits associated to cluster
  void add_cluster_size( TrackEvaluationContainerv1::ClusterStruct& cluster, TrkrCluster* trk_clus)
  {

    TrkrClusterv4 *trk_clusv4 = dynamic_cast<TrkrClusterv4*> (trk_clus);
    cluster.size = trk_clusv4->getSize();
    cluster.phi_size = trk_clusv4->getPhiSize();
    cluster.z_size = trk_clusv4->getZSize();
    cluster.ovlp = trk_clusv4->getOverlap();
    cluster.edge = trk_clusv4->getEdge();
    cluster.adc = trk_clusv4->getAdc();

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

  // ad}d truth information
  void add_truth_information( TrackEvaluationContainerv1::TrackStruct& track, PHG4Particle* particle, PHG4TruthInfoContainer* truthinfo )
  {
    if( particle )
    {
      PHG4VtxPoint* vtx  = truthinfo->GetVtx(particle->get_vtx_id());
      track.is_primary = is_primary( particle );
      track.pid = particle->get_pid();
      track.gtrackID = particle->get_track_id();
      track.truth_t = vtx->get_t();
      track.truth_px = particle->get_px();
      track.truth_py = particle->get_py();
      track.truth_pz = particle->get_pz();
      track.truth_pt = get_pt( track.truth_px, track.truth_py );
      track.truth_p = get_p( track.truth_px, track.truth_py, track.truth_pz );
      track.truth_eta = get_eta( track.truth_p, track.truth_pz );
    }
  }

  // calculate intersection between line and circle
  double line_circle_intersection( const TVector3& p0, const TVector3& p1, double radius )
  {
    const double A = square(p1.x() - p0.x()) + square(p1.y() - p0.y());
    const double B = 2*p0.x()*(p1.x()-p0.x()) + 2*p0.y()*(p1.y()-p0.y());
    const double C = square(p0.x()) + square(p0.y()) - square(radius);
    const double delta = square(B)-4*A*C;
    if( delta < 0 ) return -1;

    // check first intersection
    const double tup = (-B + std::sqrt(delta))/(2*A);
    if( tup >= 0 && tup < 1 ) return tup;

    // check second intersection
    const double tdn = (-B-sqrt(delta))/(2*A);
    if( tdn >= 0 && tdn < 1 ) return tdn;

    // no valid extrapolation
    return -1;
  }

  // print to stream
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TrackEvaluationContainerv1::ClusterStruct& cluster )
  {
    out << "ClusterStruct" << std::endl;
    out << "  cluster: (" << cluster.x << "," << cluster.y << "," << cluster.z << ")" << std::endl;
    out << "  track:   (" << cluster.trk_x << "," << cluster.trk_y << "," << cluster.trk_z << ")" << std::endl;
    out << "  truth:   (" << cluster.truth_x << "," << cluster.truth_y << "," << cluster.truth_z << ")" << std::endl;
    return out;
  }

  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TVector3& position)
  {
    out << "(" << position.x() << ", " << position.y() << ", " << position.z() << ")";
    return out;
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

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if( !m_cluster_map )
  { m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); }

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

  // tpc geometry
  m_tpc_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert( m_tpc_geom_container );

  // micromegas geometry
  m_micromegas_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void TrackEvaluation::evaluate_event()
{
  if(!m_container) return;

  // create event struct
  TrackEvaluationContainerv1::EventStruct event;
  if(m_cluster_map)
  {
    for(const auto& hitsetkey:m_cluster_map->getHitSetKeys())
    {
      const auto trkrId = TrkrDefs::getTrkrId(hitsetkey);
      const auto layer = TrkrDefs::getLayer(hitsetkey);
      assert(layer<TrackEvaluationContainerv1::EventStruct::max_layer);

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

  // store
  m_container->addEvent(event);
}

//_____________________________________________________________________
void TrackEvaluation::evaluate_clusters()
{

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  m_container->clearClusters();
  SvtxTrack *track = nullptr;
  // first loop over hitsets
  for( const auto& hitsetkey:m_cluster_map->getHitSetKeys())
  {
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // create cluster structure
      auto cluster_struct = create_cluster( key, cluster, track );
      add_cluster_size( cluster_struct, cluster);
      add_cluster_energy( cluster_struct, key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( key );
      const bool is_micromegas( TrkrDefs::getTrkrId(key) == TrkrDefs::micromegasId );
      if( is_micromegas )
      {
        const int tileid = MicromegasDefs::getTileId(key);
        add_truth_information_micromegas( cluster_struct, tileid, g4hits );
      } else {
        add_truth_information( cluster_struct, g4hits );
      }

      // add in array
      m_container->addCluster( cluster_struct );
    }
  }

}

//_____________________________________________________________________
void TrackEvaluation::evaluate_tracks()
{

  if( !( m_track_map && m_cluster_map && m_container ) ) 
  { return; }
      
  // clear array
  m_container->clearTracks();

  for( const auto& [track_id,track]:*m_track_map )
  {
    
    auto track_struct = create_track( track );
    
    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct.contributors = contributors;

    // get particle
    auto particle = m_g4truthinfo->GetParticle(id);
    track_struct.embed = get_embed(particle);
    ::add_truth_information(track_struct, particle,m_g4truthinfo);

    // running iterator over track states, used to match a given cluster to a track state
    auto state_iter = track->begin_states();

    // loop over clusters
    for( const auto& cluster_key:get_cluster_keys( track ) )
    {
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackEvaluation::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      // create new cluster struct
      auto cluster_struct = create_cluster( cluster_key, cluster, track );
      add_cluster_size( cluster_struct, cluster);
      add_cluster_energy( cluster_struct, cluster_key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( cluster_key );
      const bool is_micromegas( TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::micromegasId );
      if( is_micromegas )
      {
        const int tileid = MicromegasDefs::getTileId(cluster_key);
        add_truth_information_micromegas( cluster_struct, tileid, g4hits );
      } else {
        add_truth_information( cluster_struct, g4hits );
      }

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
      if( is_micromegas )
      {
        const int tileid = MicromegasDefs::getTileId(cluster_key);
        add_trk_information_micromegas( cluster_struct, tileid, state_iter->second );
      } else {
        add_trk_information( cluster_struct, state_iter->second );
      }

      // add to track
      track_struct.clusters.push_back( cluster_struct );
    }
    m_container->addTrack( track_struct );
  }
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
  for( const auto& cluster_key:get_cluster_keys( track ) )
  {
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

//_____________________________________________________________________
TrackEvaluationContainerv1::ClusterStruct TrackEvaluation::create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster, SvtxTrack *track ) const
{

  TrackSeed *si_seed = nullptr;
  TrackSeed *tpc_seed = nullptr;
  if(track!=nullptr){
    si_seed = track->get_silicon_seed();
    tpc_seed = track->get_tpc_seed();
  }

  // get global coordinates
  Acts::Vector3 global;
  global = m_tGeometry->getGlobalPosition(key, cluster);
  TrackEvaluationContainerv1::ClusterStruct cluster_struct;
  cluster_struct.layer = TrkrDefs::getLayer(key);
  cluster_struct.x = global.x();
  cluster_struct.y = global.y();
  cluster_struct.z = global.z();
  cluster_struct.r = get_r( cluster_struct.x, cluster_struct.y );
  cluster_struct.phi = std::atan2( cluster_struct.y, cluster_struct.x );
  cluster_struct.phi_error = 0.0;
  cluster_struct.z_error = 0.0;
  cluster_struct.trk_alpha = 0.0;
  cluster_struct.trk_beta = 0.0;
  ClusterErrorPara ClusErrPara;
  if(track!=0){
    float r = cluster_struct.r;
    if(cluster_struct.layer>7){
      auto para_errors_mm = ClusErrPara.get_cluster_error(tpc_seed,cluster,r,key);
      cluster_struct.phi_error = sqrt(para_errors_mm.first)/cluster_struct.r;
      cluster_struct.z_error = sqrt(para_errors_mm.second);
      //	float R = TMath::Abs(1.0/tpc_seed->get_qOverR());
      cluster_struct.trk_radius = 1.0/tpc_seed->get_qOverR();
      cluster_struct.trk_alpha = (r*r) /(2*r*TMath::Abs(1.0/tpc_seed->get_qOverR()));
      cluster_struct.trk_beta = atan(tpc_seed->get_slope());
    }else{
      auto para_errors_mvtx = ClusErrPara.get_cluster_error(si_seed,cluster,r,key);
      cluster_struct.phi_error = sqrt(para_errors_mvtx.first)/cluster_struct.r;
      cluster_struct.z_error = sqrt(para_errors_mvtx.second);	
      //	float R = TMath::Abs(1.0/si_seed->get_qOverR());
      cluster_struct.trk_radius = 1.0/tpc_seed->get_qOverR();
      cluster_struct.trk_alpha = (r*r) /(2*r*TMath::Abs(1.0/tpc_seed->get_qOverR()));
      cluster_struct.trk_beta = atan(si_seed->get_slope());
    }
  }
  return cluster_struct;
}

//_____________________________________________________________________
void TrackEvaluation::add_trk_information( TrackEvaluationContainerv1::ClusterStruct& cluster, SvtxTrackState* state ) const
{
  // need to extrapolate to the right r
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

//_____________________________________________________________________
void TrackEvaluation::add_trk_information_micromegas( TrackEvaluationContainerv1::ClusterStruct& cluster, int tileid, SvtxTrackState* state ) const
{

  // get geometry cylinder from layer
  const auto layer = cluster.layer;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geom_container->GetLayerGeom(layer));
  assert( layergeom );

  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster.x, cluster.y, cluster.z );
  const auto cluster_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, cluster_world );
  
  // convert track position to local tile coordinates
  TVector3 track_world( state->get_x(), state->get_y(), state->get_z() );
  TVector3 track_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, track_world );

  // convert direction to local tile coordinates
  const TVector3 direction_world( state->get_px(), state->get_py(), state->get_pz() );
  const TVector3 direction_local = layergeom->get_local_from_world_vect( tileid, m_tGeometry, direction_world );

  // extrapolate to same local z (should be zero) as cluster
  const auto delta_z = cluster_local.z() - track_local.z();
  track_local += TVector3(
    delta_z*direction_local.x()/direction_local.z(),
    delta_z*direction_local.y()/direction_local.z(),
    delta_z );

  // convert back to global coordinates
  track_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, track_local );

  // store state position
  cluster.trk_x = track_world.x();
  cluster.trk_y = track_world.y();
  cluster.trk_z = track_world.z();
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

//_____________________________________________________________________
void TrackEvaluation::add_truth_information( TrackEvaluationContainerv1::ClusterStruct& cluster, std::set<PHG4Hit*> g4hits ) const
{
  // store number of contributing g4hits
  cluster.truth_size = g4hits.size();

  // get layer, tpc flag and corresponding layer geometry
  const auto layer = cluster.layer;
  const bool is_tpc( layer >= 7 && layer < 55 );
  const PHG4TpcCylinderGeom* layergeom = is_tpc ? m_tpc_geom_container->GetLayerCellGeom(layer):nullptr;
  const auto rin = layergeom ? layergeom->get_radius()-layergeom->get_thickness()/2:0;
  const auto rout = layergeom ? layergeom->get_radius()+layergeom->get_thickness()/2:0;

  // convert hits to list of interpolation_data_t
  interpolation_data_t::list hits;
  for( const auto& g4hit:g4hits )
  {
    interpolation_data_t::list tmp_hits;
    const auto weight = g4hit->get_edep();
    for( int i = 0; i < 2; ++i )
    {
      const TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
      const TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
      tmp_hits.push_back( {.position = g4hit_world, .momentum = momentum_world, .weight = weight } );
    }

    if( is_tpc )
    {
      // add layer boundary checks
      // ensure first hit has lowest r
      auto r0 = get_r(tmp_hits[0].x(),tmp_hits[0].y());
      auto r1 = get_r(tmp_hits[1].x(),tmp_hits[1].y());
      if( r0 > r1 )
      {
        std::swap(tmp_hits[0],tmp_hits[1]);
        std::swap(r0, r1);
      }

      // do nothing if out of bound
      if( r1 <= rin || r0 >= rout ) continue;

      // keep track of original deltar
      const auto dr_old = r1-r0;

      // clamp r0 to rin
      if( r0 < rin )
      {
        const auto t = line_circle_intersection( tmp_hits[0].position, tmp_hits[1].position, rin );
        if( t<0 ) continue;

        tmp_hits[0].position = tmp_hits[0].position*(1.-t) + tmp_hits[1].position*t;
        tmp_hits[0].momentum = tmp_hits[0].momentum*(1.-t) + tmp_hits[1].momentum*t;
        r0 = rin;
      }

      if( r1 > rout )
      {
        const auto t = line_circle_intersection( tmp_hits[0].position, tmp_hits[1].position, rout );
        if( t<0 ) continue;

        tmp_hits[1].position = tmp_hits[0].position*(1.-t) + tmp_hits[1].position*t;
        tmp_hits[1].momentum = tmp_hits[0].momentum*(1.-t) + tmp_hits[1].momentum*t;
        r1 = rout;
      }

      // update weights, only if clamping occured
      const auto dr_new = r1-r0;
      tmp_hits[0].weight *= dr_new/dr_old;
      tmp_hits[1].weight *= dr_new/dr_old;
    }

    // store in global list
    hits.push_back(std::move(tmp_hits[0]));
    hits.push_back(std::move(tmp_hits[1]));
  }

  // add truth position
  cluster.truth_x = average<&interpolation_data_t::x>( hits );
  cluster.truth_y = average<&interpolation_data_t::y>( hits );
  cluster.truth_z = average<&interpolation_data_t::z>( hits );

  // add truth momentum information
  cluster.truth_px = average<&interpolation_data_t::px>( hits );
  cluster.truth_py = average<&interpolation_data_t::py>( hits );
  cluster.truth_pz = average<&interpolation_data_t::pz>( hits );

  cluster.truth_r = get_r( cluster.truth_x, cluster.truth_y );
  cluster.truth_phi = std::atan2( cluster.truth_y, cluster.truth_x );

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

}

//_____________________________________________________________________
void TrackEvaluation::add_truth_information_micromegas( TrackEvaluationContainerv1::ClusterStruct& cluster, int tileid, std::set<PHG4Hit*> g4hits ) const
{
  // store number of contributing g4hits
  cluster.truth_size = g4hits.size();

  const auto layer = cluster.layer;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geom_container->GetLayerGeom(layer));
  assert( layergeom );

  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster.x, cluster.y, cluster.z );
  const auto cluster_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, cluster_world );

  // convert hits to list of interpolation_data_t
  interpolation_data_t::list hits;
  for( const auto& g4hit:g4hits )
  {
    const auto weight = g4hit->get_edep();
    for( int i = 0; i < 2; ++i )
    {

      // convert position to local
      TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
      auto g4hit_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, g4hit_world );

      // convert momentum to local
      TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
      TVector3 momentum_local = layergeom->get_local_from_world_vect( tileid, m_tGeometry, momentum_world );

      hits.push_back( {.position = g4hit_local, .momentum = momentum_local, .weight = weight } );
    }
  }

  // do position interpolation
  const TVector3 interpolation_local(
    average<&interpolation_data_t::x>(hits),
    average<&interpolation_data_t::y>(hits),
    average<&interpolation_data_t::z>(hits) );

  // convert back to global
  const TVector3 interpolation_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, interpolation_local );

  // do momentum interpolation
  const TVector3 momentum_local(
    average<&interpolation_data_t::px>(hits),
    average<&interpolation_data_t::py>(hits),
    average<&interpolation_data_t::pz>(hits));

  // convert back to global
  const TVector3 momentum_world = layergeom->get_world_from_local_vect( tileid, m_tGeometry, momentum_local );

  cluster.truth_x = interpolation_world.x();
  cluster.truth_y = interpolation_world.y();
  cluster.truth_z = interpolation_world.z();
  cluster.truth_r = get_r( cluster.truth_x, cluster.truth_y );
  cluster.truth_phi = std::atan2( cluster.truth_y, cluster.truth_x );

  /* add truth momentum information */
  cluster.truth_px = momentum_world.x();
  cluster.truth_py = momentum_world.y();
  cluster.truth_pz = momentum_world.z();

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

}

/*!
 * \file DSTEmulation.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "DSTEmulator.h"
#include "TrackEvaluationContainerv1.h"
#include "DSTCompressor.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/recoConsts.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterv2.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <TFile.h>
#include <TNtuple.h>

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

  /*  //! radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }
  */
  /*
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
  */
  /*
  //! true if a track is a primary
  inline int is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }
  */
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
#ifdef JUNK
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
    std::cout << " (x|y|z|r|l) " 
	 << cluster_struct.x << " | " 
	 << cluster_struct.y << " | " 
	 << cluster_struct.z << " | " 
	 << cluster_struct.r << " | " 
	 << cluster_struct.layer << " | " 
	      << std::endl;
    std::cout << " (xl|yl) " 
	 << cluster->getLocalX() << " | " 
	 << cluster->getLocalY() 
	      << std::endl;
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
    std::cout << "inter" << "\n";

    cluster.truth_x = interpolate<&PHG4Hit::get_x>( hits, rextrap );
    cluster.truth_y = interpolate<&PHG4Hit::get_y>( hits, rextrap );
    cluster.truth_z = interpolate<&PHG4Hit::get_z>( hits, rextrap );
    cluster.truth_r = get_r( cluster.truth_x, cluster.truth_y );
    cluster.truth_phi = std::atan2( cluster.truth_y, cluster.truth_x );

    /* add truth momentum information */
    cluster.truth_px = interpolate<&PHG4Hit::get_px>( hits, rextrap );
    cluster.truth_py = interpolate<&PHG4Hit::get_py>( hits, rextrap );
    cluster.truth_pz = interpolate<&PHG4Hit::get_pz>( hits, rextrap );

    std::cout << "inter2" << "\n";

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
#endif
}

//_____________________________________________________________________
DSTEmulator::DSTEmulator( const std::string& name, const std::string &filename, int inBits,
                          int inSabotage, bool compress):
  SubsysReco( name)
  , _filename(filename)
  , _tfile(nullptr)
  , nBits(inBits)
  , sabotage(inSabotage)
  , apply_compression(compress)
{}

//_____________________________________________________________________
int DSTEmulator::Init(PHCompositeNode* topNode )
{ 
  if (_tfile) delete _tfile;
  _tfile = new TFile(_filename.c_str(), "RECREATE");
 
  _dst_data = new TNtuple("dst_data", "dst data","event:seed:"
			  "c3x:c3y:c3z:c3p:t3x:t3y:t3z:"
			  "c2x:c2y:c2r:c2l:t2x:t2y:t2r:t2l:"
			  "d2x:d2y:dr:"
			  "cmp_d2x:cmp_d2y:"
			  "pt:eta:phi:charge");
  
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTEmulator::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "DSTEmulator::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }
  
  auto newNode = new PHIODataNode<PHObject>( new TrackEvaluationContainerv1, "TrackEvaluationContainer","PHObject");
  evalNode->addNode(newNode);
  
  // m_compressor = new DSTCompressor(4.08407e-02,
  //                                  7.46530e-02,
  //                                  5.14381e-05,
  //                                  2.06291e-01,
  // with new distortion setting
  // m_compressor = new DSTCompressor(-2.70072e-02,
  //                                  2.49574e-02,
  //                                  1.12803e-03,
  //                                  5.91965e-02,
  // with mininum layer 7
  m_compressor = new DSTCompressor(6.96257e-04,
                                   3.16806e-02,
                                   7.32860e-05,
                                   5.93230e-02,
                                   nBits,
                                   nBits);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTEmulator::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int DSTEmulator::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,
						 "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE
		<< "ActsTrackingGeometry not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  // cleanup output container
  if( m_container ) m_container->Reset();

  evaluate_tracks();

  // clear maps
  m_g4hit_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTEmulator::End(PHCompositeNode* )
{ 
  _tfile->cd();
  _dst_data->Write();
  _tfile->Close();
  delete _tfile;
  return Fun4AllReturnCodes::EVENT_OK; 
}

Acts::Vector3 DSTEmulator::getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  Acts::Vector3 globalpos;
  globalpos = m_tGeometry->getGlobalPosition(key, cluster);
  
  // check if TPC distortion correction are in place and apply
  // if( m_dcc ) { globalpos = m_distortionCorrection.get_corrected_position( globalpos, m_dcc ); }

  return globalpos;
}

//_____________________________________________________________________
int DSTEmulator::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMapTPCOnly");
  if(m_track_map)
    {
      std::cout << " DSTEmulator: Using TPC Only Track Map node " << std::endl;
    }
  else
    {
      std::cout << " DSTEmulator: TPC Only Track Map node not found, using default" << std::endl;
      m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    }
  // cluster map

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(m_cluster_map){
    if(m_cluster_map->size()>0)
      {
	std::cout << " DSTEmulator: Using CORRECTED_TRKR_CLUSTER node " << std::endl;
      }
    else
      {
	std::cout << " DSTEmulator: CORRECTED_TRKR_CLUSTER node not found, using TRKR_CLUSTER" << std::endl;
	m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
      }
  }else{

    std::cout << " DSTEmulator: CORRECTED_TRKR_CLUSTER node not found at all, using TRKR_CLUSTER" << std::endl;
    m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  }
  

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
void DSTEmulator::evaluate_tracks()
{

  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  int iseed = 0;
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED")){
    iseed = rc->get_IntFlag("RANDOMSEED");
  }

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
    //    add_truth_information(track_struct, particle);

    // running iterator over track states, used to match a given cluster to a track state
    auto state_iter = track->begin_states();

    // loop over clusters
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

      if( !cluster )
      {
        std::cout << "DSTEmulator::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      if(TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId) continue;

      // create new cluster struct
      // auto cluster_struct = create_cluster( cluster_key, cluster );

      // find track state that is the closest to cluster
      /* this assumes that both clusters and states are sorted along r */
      //     const auto radius( cluster_struct.r );
      const Acts::Vector3 globalpos_d = getGlobalPosition(cluster_key, cluster);
      float radius = get_r( globalpos_d[0], globalpos_d[1] );
      float clu_phi = std::atan2( globalpos_d[0], globalpos_d[1] );
      std::cout << "radius " << radius << std::endl;
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
      // Got cluster and got state
      /*
       */

      std::cout << "NEW (xg|yg) " 
		<< globalpos_d[0] << " | " 
		<<  globalpos_d[1]
		<< std::endl;
      std::cout << "NEW (xl|yl) " 
		<< cluster->getLocalX() << " | " 
		<< cluster->getLocalY() 
		<< std::endl;
      
      //state to local
      // store track state in cluster struct
      //      add_trk_information( cluster_struct, state_iter->second );
      //	void add_trk_information( TrackEvaluationContainerv1::ClusterStruct& cluster, SvtxTrackState* state )
      // need to extrapolate the track state to the right cluster radius to get proper residuals  
      const auto trk_r = get_r( state_iter->second->get_x(), state_iter->second->get_y() );
      std::cout << " trk_r  " << trk_r << std::endl;
      const auto dr = get_r( globalpos_d[0], globalpos_d[1] ) - trk_r;
      std::cout << " dr  " << dr << std::endl;
      const auto trk_drdt = (state_iter->second->get_x()*state_iter->second->get_px() + state_iter->second->get_y()*state_iter->second->get_py())/trk_r;
      std::cout << " trk_drdt  " << trk_drdt << std::endl;
      const auto trk_dxdr = state_iter->second->get_px()/trk_drdt;
      std::cout << " trk_dxdr " << trk_dxdr << std::endl;
      const auto trk_dydr = state_iter->second->get_py()/trk_drdt;
      std::cout << " trk_dydr  " << trk_dydr << std::endl;
      const auto trk_dzdr = state_iter->second->get_pz()/trk_drdt;
      std::cout << " trk_dzdr  " << trk_dzdr << std::endl;
      
	// store state position
      /*	*/
      float trk_x = state_iter->second->get_x() + dr*trk_dxdr;
      float trk_y = state_iter->second->get_y() + dr*trk_dydr;
      float trk_z = state_iter->second->get_z() + dr*trk_dzdr;
      std::cout << "trk_x " << state_iter->second->get_x() << "trk_y" << state_iter->second->get_y() << "trk_z " << state_iter->second->get_z() << std::endl;
      //      float trk_r = get_r( trk_x, trk_y );
      //cluster.trk_phi = std::atan2( cluster.trk_y, cluster.trk_x );

      
      /* store local momentum information
	 cluster.trk_px = state->get_px();
	 cluster.trk_py = state->get_py();
	 cluster.trk_pz = state->get_pz();
      */
      auto layer = TrkrDefs::getLayer(cluster_key);
      /*
      std::cout << " track 3D (x|y|z|r|l) " 
		<< trk_x << " | " 
		<< trk_y << " | " 
		<< trk_z << " | " 
		<< trk_r << " | " 
		<< layer << " | " 
		<< std::endl;

       std::cout << " cluster 3D (x|y|z|r|l) " 
	      << cluster->getX() << " | " 
	      << cluster->getY() << " | " 
	      << cluster->getZ() << " | " 
	      << get_r(cluster->getX() ,cluster->getY() ) << " | " 
	      << layer << " | "
	      << std::endl;
      */
    //Get Surface
    Acts::Vector3 global(trk_x, trk_y, trk_z);
    //    TrkrDefs::subsurfkey subsurfkey = cluster->getSubSurfKey();

    //    std::cout << " subsurfkey: " << subsurfkey << std::endl;
    auto mapIter = m_tGeometry->maps().m_tpcSurfaceMap.find(layer);
    
    if(mapIter == m_tGeometry->maps().m_tpcSurfaceMap.end()){
      std::cout << PHWHERE 
		<< "Error: hitsetkey not found in clusterSurfaceMap, layer = " << 
	trk_r//layer 
		<< " hitsetkey = "
		<< hitsetkey << std::endl;
      continue;// nullptr;
    }
    std::cout << " g0: " << global[0] << " g1: " << global[1] << " g2:" << global[2] << std::endl;
    // double global_phi = atan2(global[1], global[0]);
    // double global_z = global[2];
  

    // Predict which surface index this phi and z will correspond to
    // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
    std::vector<Surface> surf_vec = mapIter->second;

    Acts::Vector3 world(globalpos_d[0], globalpos_d[1],globalpos_d[2]);
    double world_phi = atan2(world[1], world[0]);
    double world_z = world[2];

    double fraction =  (world_phi + M_PI) / (2.0 * M_PI);
    double rounded_nsurf = round( (double) (surf_vec.size()/2) * fraction  - 0.5);
    unsigned int nsurf = (unsigned int) rounded_nsurf; 
    if(world_z < 0)
      nsurf += surf_vec.size()/2;
    
    Surface surface = surf_vec[nsurf];

    Acts::Vector3 center = surface->center(m_tGeometry->geometry().getGeoContext()) 
      / Acts::UnitConstants::cm;
  
    // no conversion needed, only used in acts
    //    Acts::Vector3 normal = surface->normal(m_tGeometry->getGeoContext());
    double TrkRadius = sqrt(trk_x * trk_x + trk_y * trk_y);
    double rTrkPhi = TrkRadius * atan2(trk_y, trk_x);//trkphi;
    double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
    double surfPhiCenter = atan2(center[1], center[0]);
    double surfRphiCenter = surfPhiCenter * surfRadius;
    double surfZCenter = center[2];

    double trklocX =  0;
    double trklocY =  0;
    float delta_r = 0;

    delta_r = TrkRadius - surfRadius;
    trklocX = rTrkPhi - surfRphiCenter;
    trklocY = trk_z - surfZCenter;
    /*
    std::cout << " clus 2D: (xl|yl) " 
	      << cluster->getLocalX() << " | " 
	      << cluster->getLocalY() 
	      << std::endl;

    std::cout << " trk 2D: (xl|yl) " 
	      << trklocX << " | " 
	      << trklocY 
	      << std::endl;
    */
    float delta_x = trklocX - cluster->getLocalX();
    float delta_y = trklocY - cluster->getLocalY();
    
    float comp_dx = compress_dx(delta_x);
    float comp_dy = compress_dy(delta_y);

    // Sabotage the tracking by multiplying the residuals with a random number
    if (sabotage == -1) {
      comp_dx = rnd.Uniform(-1, 1) * delta_x;
      comp_dy = rnd.Uniform(-1, 1) * delta_y;
    } else if (sabotage == -2) {
      comp_dx = rnd.Uniform(-10, 10);
      comp_dy = rnd.Uniform(-10, 10);
    }

    /* "dst_data", "dst data","event:seed:"
			  "c3x:c3y:c3z:t3x:t3y:t3z:"
			  "c2x:c2y:c2r:c2l:t2x:t2y:t2r:t2l:"
			  "d2x:d2y:"
			  "cmp_d2x:cmp_d2y:"
			  "pt:eta:phi"
    */
    float data[] = {1,
		    (float)iseed,
		    (float)globalpos_d[0],
		    (float)globalpos_d[1],
		    (float)globalpos_d[2],
		    clu_phi,
		    trk_x,
		    trk_y,
		    trk_z,
		    cluster->getLocalX(),
		    cluster->getLocalY(),
		    get_r((float)globalpos_d[0],(float)globalpos_d[1]),
		    (float)layer,
		    (float)trklocX,
		    (float)trklocY,
		    get_r(trk_x,trk_y),
		    (float)layer,
		    delta_x,
		    delta_y,
		    delta_r,
		    comp_dx,
		    comp_dy,
		    track_struct.pt,
		    track_struct.eta,
		    track_struct.phi,
        (float) track_struct.charge
    };
    _dst_data->Fill(data);
    std::cout << "filled" << "\n";

    if (apply_compression) {
      cluster->setLocalX(trklocX - comp_dx);
      cluster->setLocalY(trklocY - comp_dy);
    }


#ifdef JUNK
    /*
    store state angles in (r,phi) and (r,z) plans
    they are needed to study space charge distortions
    */
      /*
    const auto cosphi( std::cos( cluster.trk_phi ) );
    const auto sinphi( std::sin( cluster.trk_phi ) );
    const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto trk_pz = state->get_pz();
    cluster.trk_alpha = std::atan2( trk_pphi, trk_pr );
    cluster.trk_beta = std::atan2( trk_pz, trk_pr );
    cluster.trk_phi_error = state->get_phi_error();
    cluster.trk_z_error = state->get_z_error();
      */
#endif 
    }
  }
  
  std::cout << "DSTEmulator::evaluate_tracks - tracks: " << m_container->tracks().size() << std::endl;

  
}

//_____________________________________________________________________
float DSTEmulator::compress_dx( float in_val )
{
  unsigned short key = m_compressor->compressPhi(in_val);
  float out_val = m_compressor->decompressPhi(key);
  return out_val;
}

//_____________________________________________________________________
float DSTEmulator::compress_dy( float in_val )
{
  unsigned short key = m_compressor->compressZ(in_val);
  float out_val = m_compressor->decompressZ(key);
  return out_val;
}

//_____________________________________________________________________
DSTEmulator::G4HitSet DSTEmulator::find_g4hits( TrkrDefs::cluskey cluster_key ) const
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
      else std::cout << "DSTEmulator::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> DSTEmulator::get_max_contributor( SvtxTrack* track ) const
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
int DSTEmulator::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

#include "PHTpcDeltaZCorrection.h"

#include "PHTpcDeltaZCorrection.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/TrackSeed.h>    
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <gsl/gsl_const_mksa.h> // for the speed of light
#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

namespace
{
  /// square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  /// speed of light, in cm per ns
  static constexpr double speed_of_light = GSL_CONST_MKSA_SPEED_OF_LIGHT*1e-7;

}

//____________________________________________________________________________..
PHTpcDeltaZCorrection::PHTpcDeltaZCorrection(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{ InitializeParameters(); }


//____________________________________________________________________________..
int PHTpcDeltaZCorrection::InitRun(PHCompositeNode*)
{
  UpdateParametersWithMacro();
  m_bz_const = get_double_param("bz_const");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcDeltaZCorrection::process_event(PHCompositeNode *topNode )
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHTpcDeltaZCorrection::End(PHCompositeNode* /*topNode*/ )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
void PHTpcDeltaZCorrection::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
  set_default_double_param("bz_const", 1.4);  // Tesla
  return;
}

//_____________________________________________________________________
int PHTpcDeltaZCorrection::load_nodes( PHCompositeNode* topNode )
{
  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // get necessary nodes
  m_track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  assert(m_track_map);

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if(m_cluster_map)
    {
      if(Verbosity() > 0) std::cout << " Using CORRECTED_TRKR_CLUSTER node " << std::endl;
    }
  else
    {
      if(Verbosity() > 0) std::cout << " CORRECTED_TRKR_CLUSTER node not found, using TRKR_CLUSTER" << std::endl;
      m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }
  assert(m_cluster_map);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void PHTpcDeltaZCorrection::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;
  for( unsigned int iter = 0; iter != m_track_map->size(); ++iter )
    {
      TrackSeed *seed = m_track_map->get(iter);
      if(!seed)
	{ continue; }
      process_track( iter, seed ); 
    }

  m_corrected_clusters.clear();
}

//_____________________________________________________________________
void PHTpcDeltaZCorrection::process_track( unsigned int key, TrackSeed* track )
{

  // keep track of the global position of previous cluster on track
  const Acts::Vector3 origin = {track->get_x(), track->get_y(), track->get_z()};

  // radius
  const double radius = fabs(1./track->get_qOverR()); // cm

  // helix center
  const double center_x = track->get_X0();
  const double center_y = track->get_Y0();

  // origin to center 2D vector
  const Acts::Vector2 orig_vect = {origin.x()-center_x, origin.y()-center_y };

  // print
  if( Verbosity() )
  {
    std::cout << "PHTpcDeltaZCorrection -"
      << " track: " << key
      << " positive: " << track->get_charge()
      << " center: " << center_x << ", " << center_y
      << " radius: " << radius
      << std::endl;
  }

  // loop over clusters. Assume they are ordered by layer
  for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter)
	{
    // store cluster key
    const auto& cluster_key = *key_iter;

    // consider TPC clusters only
    if( TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId ) continue;

    // skip if cluster was already corrected
    if( !m_corrected_clusters.insert(cluster_key).second ) continue;

    // get cluster
    auto cluster = m_cluster_map->findCluster( cluster_key );
    if(!cluster) continue;

    // get cluster global position
    const auto global = m_tGeometry->getGlobalPosition(cluster_key, cluster);

    // get delta z
    const double delta_z = global.z() - origin.z();

    // cluster position to center
    const Acts::Vector2 cluster_vect = {global.x()-center_x, global.y()-center_y};

    // delta phi
    const double delta_phi = std::atan2(
      cluster_vect.y()*orig_vect.x()-cluster_vect.x()*orig_vect.y(),
      cluster_vect.x()*orig_vect.x() + cluster_vect.y()*orig_vect.y() );

    // helical path length
    const double pathlength = std::sqrt( square( delta_z ) + square( radius*delta_phi ) );
    
    // adjust cluster position to account for particles propagation time
    /*
     * accounting for particles finite velocity results in reducing the electron drift time by pathlenght/c
     * this in turn affects the cluster z, so that it is always closer to the readout plane
     */
    const double t_correction = pathlength /speed_of_light;  
    cluster->setLocalY( cluster->getLocalY() - t_correction);

    if( Verbosity() )
      { std::cout << "PHTpcDeltaZCorrection::process_track - cluster: " << cluster_key 
		  << " path length: " << pathlength << " t correction " << t_correction << std::endl; }
	}
  
}

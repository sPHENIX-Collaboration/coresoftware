/**
 * \file TpcSpaceChargeReconstruction.cc
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeReconstruction.h"
#include "TpcSpaceChargeReconstructionHelper.h"
#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>

#include <cassert>
#include <memory>

namespace
{

  /// square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  /// calculate delta_phi between -pi and pi
  template< class T>
    inline constexpr T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  /// radius
  template<class T> T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  /// return number of clusters of a given type that belong to a tracks
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    return std::count_if( track->begin_cluster_keys(), track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }

  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2.*M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;

}

//_____________________________________________________________________
TpcSpaceChargeReconstruction::TpcSpaceChargeReconstruction( const std::string& name ):
  SubsysReco( name)
  , PHParameterInterface(name)
  , m_matrix_container( new TpcSpaceChargeMatrixContainerv1 )
{
  InitializeParameters();
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::set_grid_dimensions( int phibins, int rbins, int zbins )
{ m_matrix_container->set_grid_dimensions( phibins, rbins, zbins ); }

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::Init(PHCompositeNode* /*topNode*/ )
{
  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_clusters = 0;
  m_accepted_clusters = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::InitRun(PHCompositeNode* )
{

  // load parameters
  UpdateParametersWithMacro();
  m_max_talpha = get_double_param( "spacecharge_max_talpha" );
  m_max_drphi = get_double_param( "spacecharge_max_drphi" );
  m_max_tbeta = get_double_param( "spacecharge_max_tbeta" );
  m_max_dz = get_double_param( "spacecharge_max_dz" );

  // print
  if( Verbosity() )
  {
    std::cout
      << "TpcSpaceChargeReconstruction::InitRun\n"
      << " m_outputfile: " << m_outputfile << "\n"
      << " m_use_micromegas: " << std::boolalpha <<  m_use_micromegas << "\n"
      << " m_max_talpha: " << m_max_talpha << "\n"
      << " m_max_drphi: " << m_max_drphi << "\n"
      << " m_max_tbeta: " << m_max_tbeta << "\n"
      << " m_max_dz: " << m_max_dz << "\n"
      << std::endl;

    // also identify the matrix container
    m_matrix_container->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{

  // increment local event counter
  ++m_event;

  // clear global position cache
  m_globalPositions.clear();

  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::End(PHCompositeNode* /*topNode*/ )
{

  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }

  // print counters
  std::cout
    << "TpcSpaceChargeReconstruction::End -"
    << " track statistics total: " << m_total_tracks
    << " accepted: " << m_accepted_tracks
    << " fraction: " << 100.*m_accepted_tracks/m_total_tracks << "%"
    << std::endl;

  std::cout
    << "TpcSpaceChargeReconstruction::End -"
    << " cluster statistics total: " << m_total_clusters
    << " accepted: " << m_accepted_clusters << " fraction: "
    << 100.*m_accepted_clusters/m_total_clusters << "%"
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void TpcSpaceChargeReconstruction::SetDefaultParameters()
{
  // residual cuts
  set_default_double_param( "spacecharge_max_talpha", 0.6 );
  set_default_double_param( "spacecharge_max_drphi", 0.5 );
  set_default_double_param( "spacecharge_max_tbeta", 1.5 );
  set_default_double_param( "spacecharge_max_dz", 0.5 );
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // get necessary nodes
  m_tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  assert( m_tgeometry );

  m_surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  assert( m_surfmaps );

  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

    m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(m_cluster_map)
  {

    if( m_event < 2 )
    { std::cout << "TpcSpaceChargeReconstruction::load_nodes - Using CORRECTED_TRKR_CLUSTER node " << std::endl; }

  } else {

    if( m_event < 2 )
    { std::cout << "TpcSpaceChargeReconstruction::load_nodes - CORRECTED_TRKR_CLUSTER node not found, using TRKR_CLUSTER" << std::endl; }
    m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");

  }

  assert( m_cluster_map );
  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________
Acts::Vector3 TpcSpaceChargeReconstruction::get_global_position(TrkrDefs::cluskey key, TrkrCluster* cluster )
{

  // find closest iterator in map
  auto it = m_globalPositions.lower_bound( key );
  if (it == m_globalPositions.end()|| (key < it->first ))
  {
    // get global position from Acts transform
    const auto globalpos = m_transform.getGlobalPosition(
      key, 
      cluster,
      m_surfmaps,
      m_tgeometry);

    /*
     * todo: should also either apply distortion corrections
     * or make sure clusters are loaded from corrected map, after cluster mover
     */

    // add new cluster and set its key
    it = m_globalPositions.insert(it, std::make_pair(key, globalpos));
  }
  return it->second;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;

  for(const auto &[trackKey, track] : *m_track_map)
  {
    ++m_total_tracks;
    if( accept_track( track ) )
    {
      ++m_accepted_tracks;
      process_track( track );
    }
  }
}

//_____________________________________________________________________
bool TpcSpaceChargeReconstruction::accept_track( SvtxTrack* track ) const
{
  // ignore tracks whose transverse momentum is too small
  const auto pt = std::sqrt( square( track->get_px() ) + square( track->get_py() ) );
  if( pt < 0.5 ) return false;

  // ignore tracks with too few mvtx, intt and micromegas hits
  if( get_clusters<TrkrDefs::mvtxId>(track) < 2 ) return false;
  if( get_clusters<TrkrDefs::inttId>(track) < 2 ) return false;
  if( m_use_micromegas && get_clusters<TrkrDefs::micromegasId>(track) < 2 ) return false;

  // all tests passed
  return true;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_track( SvtxTrack* track )
{

  // printout all track state
  if( Verbosity() )
  {
    for( auto&& iter = track->begin_states(); iter != track->end_states(); ++iter )
    {
      const auto& [pathlength, state] = *iter;
      const auto r = std::sqrt( square( state->get_x() ) + square( state->get_y() ));
      const auto phi = std::atan2( state->get_y(), state->get_x() );
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
        << " pathlength: " << pathlength
        << " radius: " << r
        << " phi: " << phi
        << " z: " << state->get_z()
        << std::endl;
    }
  }

  // running track state
  auto state_iter = track->begin_states();

  // loop over clusters
  for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
  {

    const auto& cluster_key = *key_iter;
    auto cluster = m_cluster_map->findCluster( cluster_key );
    if( !cluster )
    {
      std::cout << PHWHERE << " unable to find cluster for key " << cluster_key << std::endl;
      continue;
    }

    ++m_total_clusters;

    // make sure
    const auto detId = TrkrDefs::getTrkrId(cluster_key);
    if(detId != TrkrDefs::tpcId) continue;

    // cluster r, phi and z
    const auto global_position = get_global_position( cluster_key, cluster );
    const auto cluster_r = get_r( global_position.x(), global_position.y() );
    const auto cluster_phi = std::atan2( global_position.y(), global_position.x() );
    const auto cluster_z = global_position.z();

    // cluster errors
    const auto cluster_rphi_error = cluster->getRPhiError();
    const auto cluster_z_error = cluster->getZError();

    /*
    remove clusters with too small errors since they are likely pathological
    and have a large contribution to the chisquare
    TODO: make these cuts configurable
    */
    if( cluster_rphi_error < 0.015 ) continue;
    if( cluster_z_error < 0.05 ) continue;

    // find track state that is the closest to cluster
    /* this assumes that both clusters and states are sorted along r */
    float dr_min = -1;
    for( auto iter = state_iter; iter != track->end_states(); ++iter )
    {
      const auto dr = std::abs( cluster_r - get_r( iter->second->get_x(), iter->second->get_y() ) );
      if( dr_min < 0 || dr < dr_min )
      {
        state_iter = iter;
        dr_min = dr;
      } else break;
    }


    // get relevant track state
    const auto state = state_iter->second;

    // extrapolate track parameters to the cluster r
    const auto track_r = get_r( state->get_x(), state->get_y() );
    const auto dr = cluster_r - track_r;
    const auto track_drdt = (state->get_x()*state->get_px() + state->get_y()*state->get_py())/track_r;
    const auto track_dxdr = state->get_px()/track_drdt;
    const auto track_dydr = state->get_py()/track_drdt;
    const auto track_dzdr = state->get_pz()/track_drdt;

    // store state position
    const auto track_x = state->get_x() + dr*track_dxdr;
    const auto track_y = state->get_y() + dr*track_dydr;
    const auto track_z = state->get_z() + dr*track_dzdr;
    const auto track_phi = std::atan2( track_y, track_x );

    // get track angles
    const auto cosphi( std::cos( track_phi ) );
    const auto sinphi( std::sin( track_phi ) );
    const auto track_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto track_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto track_pz = state->get_pz();
    const auto talpha = -track_pphi/track_pr;
    const auto tbeta = -track_pz/track_pr;

    // sanity check
    if( std::isnan(talpha) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - talpha is nan" << std::endl;
      continue;
    }

    if( std::isnan(tbeta) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - tbeta is nan" << std::endl;
      continue;
    }

    // check against limits
    if( std::abs( talpha ) > m_max_talpha ) continue;
    if( std::abs( tbeta ) > m_max_tbeta ) continue;

    // track errors
    const auto track_rphi_error = state->get_rphi_error();
    const auto track_z_error = state->get_z_error();

    // residuals
    const auto drp = cluster_r*delta_phi( cluster_phi - track_phi );
    const auto dz = cluster_z - track_z;

    // sanity checks
    if( std::isnan(drp) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - drp is nan" << std::endl;
      continue;
    }

    if( std::isnan(dz) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - dz is nan" << std::endl;
      continue;
    }

    // check against limits
    if( std::abs( drp ) > m_max_drphi ) continue;
    if( std::abs( dz ) > m_max_dz ) continue;

    // residual errors squared
    const auto erp = square(track_rphi_error) + square(cluster_rphi_error);
    const auto ez = square(track_z_error) + square(cluster_z_error);

    // sanity check
    // TODO: check whether this happens and fix upstream
    if( std::isnan( erp ) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - erp is nan" << std::endl;
      continue;
    }

    if( std::isnan( ez ) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - ez is nan" << std::endl;
      continue;
    }

    // get cell
    const auto i = get_cell_index( cluster_key, cluster );
    if( i < 0 )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - invalid cell index" << std::endl;
      continue;
    }

    // update matrices
    // see https://indico.bnl.gov/event/7440/contributions/43328/attachments/31334/49446/talk.pdf for details
    m_matrix_container->add_to_lhs(i, 0, 0, 1./erp );
    m_matrix_container->add_to_lhs(i, 0, 1, 0 );
    m_matrix_container->add_to_lhs(i, 0, 2, talpha/erp );

    m_matrix_container->add_to_lhs(i, 1, 0, 0 );
    m_matrix_container->add_to_lhs(i, 1, 1, 1./ez );
    m_matrix_container->add_to_lhs(i, 1, 2, tbeta/ez );

    m_matrix_container->add_to_lhs(i, 2, 0, talpha/erp );
    m_matrix_container->add_to_lhs(i, 2, 1, tbeta/ez );
    m_matrix_container->add_to_lhs(i, 2, 2, square(talpha)/erp + square(tbeta)/ez );

    m_matrix_container->add_to_rhs(i, 0, drp/erp );
    m_matrix_container->add_to_rhs(i, 1, dz/ez );
    m_matrix_container->add_to_rhs(i, 2, talpha*drp/erp + tbeta*dz/ez );

    // update entries in cell
    m_matrix_container->add_to_entries(i);

    // increment number of accepted clusters
    ++m_accepted_clusters;

  }

}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell_index(TrkrDefs::cluskey key, TrkrCluster* cluster )
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );


  // get global position
  const auto global_position = get_global_position( key, cluster );

  // phi
  // bound check
  float phi = std::atan2( global_position.y(), global_position.x() );
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  int iphi = phibins*(phi-m_phimin)/(m_phimax-m_phimin);

  // radius
  const float r = get_r( global_position.x(), global_position.y() );
  if( r < m_rmin || r >= m_rmax ) return -1;
  int ir = rbins*(r-m_rmin)/(m_rmax-m_rmin);

  // z
  const float z = global_position.z();
  if( z < m_zmin || z >= m_zmax ) return -1;
  int iz = zbins*(z-m_zmin)/(m_zmax-m_zmin);

  return m_matrix_container->get_cell_index( iphi, ir, iz );
}

#include "TpcSpaceChargeReconstruction.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>
#include <TH3.h>

#include <memory>

namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// calculate delta_phi between -pi and pi
  template< class T>
    T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  /// radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  /// return number of clusters of a given type that belong to a tracks
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    return std::count_if( track->begin_cluster_keys(), track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }


}

//_____________________________________________________________________
TpcSpaceChargeReconstruction::TpcSpaceChargeReconstruction( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::set_grid_dimensions( int phibins, int rbins, int zbins )
{
  m_phibins = phibins;
  m_rbins = rbins;
  m_zbins = zbins;
  m_totalbins = phibins*rbins*zbins;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::set_outputfile( const std::string& filename )
{ m_outputfile = filename; }

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::Init(PHCompositeNode* topNode )
{
  // resize vectors
  m_lhs = std::vector<matrix_t>( m_totalbins, matrix_t::Zero() );
  m_rhs = std::vector<column_t>( m_totalbins, column_t::Zero() );
  m_cluster_count = std::vector<int>( m_totalbins, 0 );
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::InitRun(PHCompositeNode* )
{
  std::cout
    << "TpcSpaceChargeReconstruction::InitRun\n"
    << " m_outputfile: " << m_outputfile << "\n"
    << " m_use_micromegas: " << std::boolalpha <<  m_use_micromegas << "\n"
    << " m_phibins: " << m_phibins << "\n"
    << " m_rbins: " << m_rbins << "\n"
    << " m_zbins: " << m_zbins << "\n"
    << " m_totalbins: " << m_totalbins << "\n"
    << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::End(PHCompositeNode* topNode )
{
  calculate_distortions( topNode );
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { if( accept_track( iter->second ) ) process_track( iter->second ); }
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

    // make sure
    const auto detId = TrkrDefs::getTrkrId(cluster_key);
    if(detId != TrkrDefs::tpcId) continue;

    // cluster r, phi and z
    const auto cluster_r = get_r( cluster->getX(), cluster->getY() );
    const auto cluster_phi = std::atan2( cluster->getY(), cluster->getX() );
    const auto cluster_z = cluster->getZ();

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

    // track errors
    #if 0
    const auto track_rphi_error = state->get_rphi_error();
    const auto track_z_error = state->get_z_error();

    // also cut on track errors
    // warning: smaller errors are probably needed when including outer tracker
    if( track_rphi_error < 0.015 ) continue;
    if( track_z_error < 0.1 ) continue;
    #endif

    // extrapolate track parameters to the cluster r
    const auto track_r = get_r( state->get_x(), state->get_y() );
    const auto dr = cluster_r - track_r;
    const auto track_drdt = get_r( state->get_px(), state->get_py() );
    const auto track_dxdr = state->get_px()/track_drdt;
    const auto track_dydr = state->get_py()/track_drdt;
    const auto track_dzdr = state->get_pz()/track_drdt;

    // store state position
    const auto track_x = state->get_x() + dr*track_dxdr;
    const auto track_y = state->get_y() + dr*track_dydr;
    const auto track_z = state->get_z() + dr*track_dzdr;
    const auto track_phi = std::atan2( track_y, track_x );

    // get residual errors squared
    // for now we only consider the error on the cluster. Not on the track
    //     const auto erp = square(track_rphi_error) + square(cluster_rphi_error);
    //     const auto ez = square(track_z_error) + square(cluster_z_error);
    const auto erp = square(cluster_rphi_error);
    const auto ez = square(cluster_z_error);

    // sanity check
    // TODO: check whether this happens and fix upstream
    if( std::isnan( erp ) ) continue;
    if( std::isnan( ez ) ) continue;

    // get residuals
    const auto drp = cluster_r*delta_phi( track_phi - cluster_phi );
    const auto dz = track_z - cluster_z;

    // sanity checks
    // TODO: check whether this happens and fix upstream
    if( std::isnan(drp) ) continue;
    if( std::isnan(dz) ) continue;

    // get track angles
    const auto cosphi( std::cos( track_phi ) );
    const auto sinphi( std::sin( track_phi ) );
    const auto track_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto track_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto track_pz = state->get_pz();
    const auto talpha = -track_pphi/track_pr;
    const auto tbeta = -track_pz/track_pr;

    // sanity check
    // TODO: check whether this happens and fix upstream
    if( std::isnan(talpha) ) continue;
    if( std::isnan(tbeta) ) continue;

    // check against limits
    static constexpr float max_talpha = 0.6;
    static constexpr float max_residual_drphi = 0.5;
    if( std::abs( talpha ) > max_talpha ) continue;
    if( std::abs( drp ) > max_residual_drphi ) continue;

    static constexpr float max_tbeta = 1.5;
    static constexpr float max_residual_dz = 0.5;
    if( std::abs( tbeta ) > max_tbeta ) continue;
    if( std::abs( dz ) > max_residual_dz ) continue;

    // get cell
    const auto i = get_cell( cluster );

    if( i < 0 || i >= m_totalbins ) continue;

    m_lhs[i](0,0) += 1./erp;
    m_lhs[i](0,1) += 0;
    m_lhs[i](0,2) += talpha/erp;

    m_lhs[i](1,0) += 0;
    m_lhs[i](1,1) += 1./ez;
    m_lhs[i](1,2) += tbeta/ez;

    m_lhs[i](2,0) += talpha/erp;
    m_lhs[i](2,1) += tbeta/ez;
    m_lhs[i](2,2) += square(talpha)/erp + square(tbeta)/ez;

    m_rhs[i](0,0) += drp/erp;
    m_rhs[i](1,0) += dz/ez;
    m_rhs[i](2,0) += talpha*drp/erp + tbeta*dz/ez;

    ++m_cluster_count[i];

  }

}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::calculate_distortions( PHCompositeNode* topNode )
{

  // create output histograms
  auto hentries = new TH3F( "hentries_rec", "hentries_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hphi = new TH3F( "hDistortionP_rec", "hDistortionP_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hz = new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );
  auto hr = new TH3F( "hDistortionR_rec", "hDistortionR_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax );

  // set axis labels
  for( auto h:{ hentries, hphi, hz, hr } )
  {
    h->GetXaxis()->SetTitle( "#phi (rad)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "z (cm)" );
  }

  // loop over bins
  for( int iphi = 0; iphi < m_phibins; ++iphi )
    for( int ir = 0; ir < m_rbins; ++ir )
    for( int iz = 0; iz < m_zbins; ++iz )
  {

    const auto icell = get_cell( iphi, ir, iz );

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 10;
    if( m_cluster_count[icell] < min_cluster_count ) continue;

    // calculate result using linear solving
    const auto cov = m_lhs[icell].inverse();
    const auto result = m_lhs[icell].partialPivLu().solve( m_rhs[icell] );

    // fill histograms
    hentries->SetBinContent( iphi+1, ir+1, iz+1, m_cluster_count[icell] );

    hphi->SetBinContent( iphi+1, ir+1, iz+1, result(0) );
    hphi->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(0,0) ) );

    hz->SetBinContent( iphi+1, ir+1, iz+1, result(1) );
    hz->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(1,1) ) );

    hr->SetBinContent( iphi+1, ir+1, iz+1, result(2) );
    hr->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(2,2) ) );

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - inverting bin " << iz << ", " << ir << ", " << iphi << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - entries: " << m_cluster_count[icell] << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - lhs: " << m_lhs[icell] << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - rhs: " << m_rhs[icell] << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - drphi: " << result(0) << " +/- " << std::sqrt( cov(0,0) ) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - dz: " << result(1) << " +/- " << std::sqrt( cov(1,1) ) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - dr: " << result(2) << " +/- " << std::sqrt( cov(2,2) ) << std::endl;
      std::cout << std::endl;
    }
  }

  // save to root file
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    for( TObject* value:{ hentries, hphi, hz, hr } ) { value->Write(); }
    outputfile->Close();
  }
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell( int iphi, int ir, int iz ) const
{
  if( iphi < 0 || iphi >= m_phibins ) return -1;
  if( ir < 0 || ir >= m_rbins ) return -1;
  if( iz < 0 || iz >= m_zbins ) return -1;
  return iz + m_zbins*( ir + m_rbins*iphi );
}

//_________________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell( float phi, float r, float z ) const
{

  // phi
  // bound check
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  int iphi = m_phibins*(phi-m_phimin)/(m_phimax-m_phimin);

  // radius
  if( r < m_rmin || r >= m_rmax ) return -1;
  int ir = m_rbins*(r-m_rmin)/(m_rmax-m_rmin);

  // z
  if( z < m_zmin || z >= m_zmax ) return -1;
  int iz = m_zbins*(z-m_zmin)/(m_zmax-m_zmin);

  return get_cell( iphi, ir, iz );
}


//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell( TrkrCluster* cluster ) const
{
  // get cluster radial coordinates
  const auto phi = std::atan2( cluster->getY(), cluster->getX() );
  const auto r = get_r( cluster->getX(), cluster->getY() );
  const auto z = cluster->getZ();
  return get_cell( phi, r, z );
}

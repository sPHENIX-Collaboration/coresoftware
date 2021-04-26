#include "TpcSpaceChargeReconstruction.h"
#include "TpcSpaceChargeReconstructionHelper.h"

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
  , PHParameterInterface(name)
{ InitializeParameters(); }

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::set_grid_dimensions( int phibins, int rbins, int zbins )
{
  m_phibins = phibins;
  m_rbins = rbins;
  m_zbins = zbins;
  m_totalbins = m_phibins*m_rbins*m_zbins;
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
  std::cout
    << "TpcSpaceChargeReconstruction::InitRun\n"
    << " m_outputfile: " << m_outputfile << "\n"
    << " m_use_micromegas: " << std::boolalpha <<  m_use_micromegas << "\n"
    << " m_phibins: " << m_phibins << "\n"
    << " m_rbins: " << m_rbins << "\n"
    << " m_zbins: " << m_zbins << "\n"
    << " m_totalbins: " << m_totalbins << "\n"
    << " m_max_talpha: " << m_max_talpha << "\n"
    << " m_max_drphi: " << m_max_drphi << "\n"
    << " m_max_tbeta: " << m_max_tbeta << "\n"
    << " m_max_dz: " << m_max_dz << "\n"
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
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert(m_cluster_map);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  {
    ++m_total_tracks;
    if( accept_track( iter->second ) )
    {
      ++m_accepted_tracks;
      process_track( iter->second );
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
      std::cout << "TpcSpaceChargeReconstruction::process_track - talpha nan" << std::endl;
      continue;
    }

    if( std::isnan(tbeta) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - tbeta nan" << std::endl;
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
      std::cout << "TpcSpaceChargeReconstruction::process_track - drp nan" << std::endl;
      continue;
    }

    if( std::isnan(dz) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - dz nan" << std::endl;
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
      std::cout << "TpcSpaceChargeReconstruction::process_track - erp nan" << std::endl;
      continue;
    }

    if( std::isnan( ez ) )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - ez nan" << std::endl;
      continue;
    }

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

    ++m_accepted_clusters;
    ++m_cluster_count[i];

  }

}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::calculate_distortions( PHCompositeNode* topNode )
{

  // create output histograms
  auto hentries( new TH3F( "hentries_rec", "hentries_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax ) );
  auto hphi( new TH3F( "hDistortionP_rec", "hDistortionP_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax ) );
  auto hz( new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax ) );
  auto hr( new TH3F( "hDistortionR_rec", "hDistortionR_rec", m_phibins, m_phimin, m_phimax, m_rbins, m_rmin, m_rmax, m_zbins, m_zmin, m_zmax ) );

  // set axis labels
  for( const auto& h:{ hentries, hphi, hz, hr } )
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

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - inverting bin " << iz << ", " << ir << ", " << iphi << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - entries: " << m_cluster_count[icell] << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - lhs: \n" << m_lhs[icell] << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - rhs: \n" << m_rhs[icell] << std::endl;
    }

    // calculate result using linear solving
    const auto cov = m_lhs[icell].inverse();
    auto partialLu = m_lhs[icell].partialPivLu();
    const auto result = partialLu.solve( m_rhs[icell] );

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
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - drphi: " << result(0) << " +/- " << std::sqrt( cov(0,0) ) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - dz: " << result(1) << " +/- " << std::sqrt( cov(1,1) ) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::calculate_distortions - dr: " << result(2) << " +/- " << std::sqrt( cov(2,2) ) << std::endl;
      std::cout << std::endl;
    }
  }

  // save everything to root file
  std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
  outputfile->cd();

  // when using migromegas, one needs to extrapolate to the rest of the acceptance
  if( m_use_micromegas )
  {
    for( const auto& h: {hentries, hphi, hr, hz} )
    {
      if( !h ) continue;
      TpcSpaceChargeReconstructionHelper::extrapolate_z(h);
      TpcSpaceChargeReconstructionHelper::extrapolate_phi1(h);
      TpcSpaceChargeReconstructionHelper::extrapolate_phi2(h);
    }
  }

  // write source histograms
  for( const auto& h: { hentries, hphi, hr, hz } ) { h->Write(); }
  
  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = []( TH3* h, const TString& name )
  {
    TH3* hneg;
    TH3* hpos;
    std::tie( hneg, hpos ) = TpcSpaceChargeReconstructionHelper::split( h );
    hneg->Write();
    hpos->Write();
    TpcSpaceChargeReconstructionHelper::copy_histogram( h, name )->Write();
    TpcSpaceChargeReconstructionHelper::copy_histogram( hneg, Form( "%s_negz", name.Data() ) )->Write();
    TpcSpaceChargeReconstructionHelper::copy_histogram( hpos, Form( "%s_posz", name.Data() ) )->Write();
  };
  
  process_histogram( hentries, "hentries" );
  process_histogram( hphi, "hIntDistortionP" );
  process_histogram( hr, "hIntDistortionR" );
  process_histogram( hz, "hIntDistortionZ" );

  // close output file
  outputfile->Close();

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

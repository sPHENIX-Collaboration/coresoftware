/**
 * \file TpcSpaceChargeReconstruction.cc
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeReconstruction.h"
#include "TpcSpaceChargeReconstructionHelper.h"
#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

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
  
  /// get sector median angle associated to a given index
  /** this assumes that sector 0 is centered on phi=0, then numbered along increasing phi */
  inline constexpr double get_sector_phi( int isec ) 
  { return isec*M_PI/6; }

  // specify bins for which one will save histograms
  static const std::vector<float> phi_rec = { get_sector_phi(9) };
  static const std::vector<float> z_rec = { 5. };

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

  /// get cluster keys from a given track
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

  /// return number of clusters of a given type that belong to a tracks
  template<int type>
    int count_clusters( const std::vector<TrkrDefs::cluskey>& keys )
  {
    return std::count_if( keys.begin(), keys.end(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }

  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const Acts::Vector3& vector )
  { 
    out << "(" << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }
  
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

  // histogram evaluation
  if( m_savehistograms ) create_histograms();

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
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_outputfile: " << m_outputfile << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_use_micromegas: " << std::boolalpha <<  m_use_micromegas << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_talpha: " << m_max_talpha << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_drphi: " << m_max_drphi << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_tbeta: " << m_max_tbeta << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_dz: " << m_max_dz << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_min_pt: " << m_min_pt << " GeV/c" << std::endl;

  // also identify the matrix container
  m_matrix_container->identify();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{

  // increment local event counter
  ++m_event;

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
  
  // save histograms
  if( m_savehistograms && m_histogramfile )
  {
    m_histogramfile->cd();    
    for( const auto& [cell,h]:m_h_drphi ) { if(h) h->Write(); }
    for( const auto& [cell,h]:m_h_dz ) { if(h) h->Write(); }
    for( const auto& [cell,h]:m_h_drphi_alpha ) { if(h) h->Write(); }
    for( const auto& [cell,h]:m_h_dz_beta ) { if(h) h->Write(); }
    m_histogramfile->Close();
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
  m_tgeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  assert( m_tgeometry );

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
  
  
  // tpc distortion corrections
  m_dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  m_dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  m_dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::create_histograms()
{
  std::cout << "TpcSpaceChargeReconstruction::create_histograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );
  
  // get bins corresponding to selected angles
  std::set<int> phibin_rec;
  std::transform( phi_rec.begin(), phi_rec.end(), std::inserter( phibin_rec, phibin_rec.end() ), [&]( const float& phi ) { return phibins*(phi-m_phimin)/(m_phimax-m_phimin); } );
  
  std::set<int> zbin_rec;
  std::transform( z_rec.begin(), z_rec.end(), std::inserter( zbin_rec, zbin_rec.end() ), [&]( const float& z ) { return zbins*(z-m_zmin)/(m_zmax-m_zmin); } );
  
  // keep track of all cell ids that match selected histograms
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {
    
    if( phibin_rec.find( iphi ) == phibin_rec.end() || zbin_rec.find( iz ) == zbin_rec.end() ) continue;
    const auto icell = m_matrix_container->get_cell_index( iphi, ir, iz );
    
    {
      // rphi residuals
      const auto hname = Form( "residual_drphi_p%i_r%i_z%i", iphi, ir, iz );
      auto h = new TH1F( hname, hname, 100, -m_max_drphi, +m_max_drphi );
      h->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
      m_h_drphi.insert( std::make_pair( icell, h ) );
    }
    
    {
      // 2D histograms
      const auto hname = Form( "residual_2d_drphi_p%i_r%i_z%i", iphi, ir, iz );
      auto h = new TH2F( hname, hname, 100, -m_max_talpha, m_max_talpha, 100, -m_max_drphi, +m_max_drphi );
      h->GetXaxis()->SetTitle( "tan#alpha" );
      h->GetYaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
      m_h_drphi_alpha.insert( std::make_pair( icell, h ) );
    }
    
    {
      // z residuals
      const auto hname = Form( "residual_dz_p%i_r%i_z%i", iphi, ir, iz );
      auto h = new TH1F( hname, hname, 100, -m_max_dz, +m_max_dz );
      h->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
      m_h_dz.insert( std::make_pair( icell, h ) );
    }
    
    {
      // 2D histograms
      static constexpr double max_tbeta = 0.5;
      const auto hname = Form( "residual_2d_dz_p%i_r%i_z%i", iphi, ir, iz );
      auto h = new TH2F( hname, hname, 100, -max_tbeta, max_tbeta, 100, -m_max_dz, +m_max_dz );
      h->GetXaxis()->SetTitle( "tan#beta" );
      h->GetYaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
      m_h_dz_beta.insert( std::make_pair( icell, h ) );
    }     
  }  
}

//_________________________________________________________________________________
Acts::Vector3 TpcSpaceChargeReconstruction::get_global_position(TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing ) const
{
  // get global position from Acts transform
  auto globalPosition = m_tgeometry->getGlobalPosition(key, cluster);
  
  // for the TPC calculate the proper z based on crossing and side
  const auto trkrid = TrkrDefs::getTrkrId(key);
  if(trkrid ==  TrkrDefs::tpcId)
  {	 
    const auto side = TpcDefs::getSide(key);
    globalPosition.z() = m_clusterCrossingCorrection.correctZ(globalPosition.z(), side, crossing);    
        
    // apply distortion corrections
    if(m_dcc_static) 
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_static ); 
    }
    
    if(m_dcc_average) 
    { 
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_average ); 
    }
    
    if(m_dcc_fluctuation) 
    { 
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_fluctuation ); 
    }
  }
    
  return globalPosition;
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

  // track pt
  if(track->get_pt() < m_min_pt)
  { return false; }

  // ignore tracks with too few mvtx, intt and micromegas hits
  const auto cluster_keys( get_cluster_keys( track ) );
  if( count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 2 ) return false;
  if( count_clusters<TrkrDefs::inttId>(cluster_keys) < 2 ) return false;
  if( m_use_micromegas && count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2 ) return false;

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

  // store crossing
  // it is needed for geting cluster's global position
  const auto crossing = track->get_crossing();
  assert( crossing != SHRT_MAX );

  // running track state
  auto state_iter = track->begin_states();
  
  // loop over clusters
  for( const auto& cluster_key:get_cluster_keys( track ) )
  {

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
    const auto global_position = get_global_position( cluster_key, cluster, crossing );
    const auto cluster_r = get_r( global_position.x(), global_position.y() );
    const auto cluster_phi = std::atan2( global_position.y(), global_position.x() );
    const auto cluster_z = global_position.z();

    // cluster errors
    double cluster_rphi_error = 0;
    double cluster_z_error = 0;
    if( m_cluster_version >= 4 )
    {
      const auto errors_square = m_cluster_error_parametrization.get_cluster_error( track->get_tpc_seed(), cluster, cluster_r, cluster_key ); 
      cluster_rphi_error = std::sqrt( errors_square.first );
      cluster_z_error = std::sqrt( errors_square.second );
    } else {
      cluster_rphi_error = cluster->getRPhiError();
      cluster_z_error = cluster->getZError();
    }

    /* 
     * as instructed by Christof, it should not be necessary to cut on small
     * cluster errors any more with clusters of version 4 or higher 
     */ 
    if( m_cluster_version < 4 )
    {
      /*
       * remove clusters with too small errors since they are likely pathological
       * and have a large contribution to the chisquare
       * TODO: make these cuts configurable
       */
      if( cluster_rphi_error < 0.015 ) continue;
      if( cluster_z_error < 0.05 ) continue;
    }
    
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
    const auto i = get_cell_index( global_position );
    if( i < 0 )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - invalid cell index" << std::endl;
      continue;
    }
    
    if( m_savehistograms )
    {      
      { const auto iter = m_h_drphi.find( i ); if( iter != m_h_drphi.end() ) iter->second->Fill( drp ); }
      { const auto iter = m_h_drphi_alpha.find( i ); if( iter != m_h_drphi_alpha.end() ) iter->second->Fill( talpha, drp ); }
      { const auto iter = m_h_dz.find( i ); if( iter != m_h_dz.end() ) iter->second->Fill( dz ); }
      { const auto iter = m_h_dz_beta.find( i ); if( iter != m_h_dz_beta.end() ) iter->second->Fill( tbeta, dz ); }
    }
    
    // check against limits
    if( std::abs( talpha ) > m_max_talpha ) continue;
    if( std::abs( tbeta ) > m_max_tbeta ) continue;

    // check against limits
    if( std::abs( drp ) > m_max_drphi ) continue;
    if( std::abs( dz ) > m_max_dz ) continue;

    if( Verbosity() )
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - layer: " << (int) TrkrDefs::getLayer(cluster_key) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
        << " cluster: (" << cluster_r << ", " << cluster_r*cluster_phi << ", " << cluster_z << ")"
        << " (" << cluster_rphi_error << ", " << cluster_z_error << ")" 
        << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
        << " track: (" << track_r << ", " << cluster_r*track_phi << ", " << track_z << ")"
        << " (" << talpha << ", " << tbeta << ")"
        << " (" << track_rphi_error << ", " << track_z_error << ")"
        << std::endl;
      std::cout << std::endl;
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
int TpcSpaceChargeReconstruction::get_cell_index( const Acts::Vector3& global_position )
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );

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

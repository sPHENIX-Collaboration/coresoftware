#include "PHSpaceChargeReconstruction.h"

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
#include <TGraphErrors.h>

#include <memory>

namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// calculate delta_phi between -pi and pi
  template< class T>
    T delta_phi( const T& phi )
  {
    if( phi > M_PI ) return phi - 2*M_PI;
    else if( phi <= -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  /// radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  /// phi
  template<class T> T get_phi( T x, T y ) { return std::atan2( y, x ); }
}

//_____________________________________________________________________
PHSpaceChargeReconstruction::PHSpaceChargeReconstruction( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::set_tpc_layers( unsigned int first_layer, unsigned int n_layers )
{
  m_firstlayer_tpc = first_layer;
  m_nlayers_tpc = n_layers;
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::set_grid_dimensions( int zbins, int rbins, int phibins )
{
  m_zbins = zbins;
  m_rbins = rbins;
  m_phibins = phibins;
  m_totalbins = zbins*rbins*phibins;
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::set_outputfile( const std::string& filename )
{ m_outputfile = filename; }

//_____________________________________________________________________
int PHSpaceChargeReconstruction::Init(PHCompositeNode* topNode )
{

  // resize vectors
  m_lhs = std::vector<matrix_t>( m_totalbins, matrix_t::Zero() );
  m_rhs = std::vector<column_t>( m_totalbins, column_t::Zero() );
  m_cluster_count = std::vector<int>( m_totalbins, 0 );
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int PHSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::End(PHCompositeNode* topNode )
{
  calculate_distortions( topNode );
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;

  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::process_track( SvtxTrack* track )
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

    // should make sure that cluster belongs to TPC
    const auto layer = TrkrDefs::getLayer(cluster_key);
    if( layer < m_firstlayer_tpc || layer >= m_firstlayer_tpc + m_nlayers_tpc ) continue;

    // cluster r, phi and z
    const auto cluster_r = get_r( cluster->getX(), cluster->getY() );
    const auto cluster_phi = get_phi( cluster->getX(), cluster->getY() );
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
    const auto track_rphi_error = state->get_rphi_error();
    const auto track_z_error = state->get_z_error();

    // also cut on track errors
    // warning: smaller errors are probably needed when including outer tracker
    if( track_rphi_error < 0.015 ) continue;
    if( track_z_error < 0.1 ) continue;

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
    const auto track_phi = get_phi( track_x, track_y );

    // get residual errors squared
    const auto erp = square(track_rphi_error) + square(cluster_rphi_error);
    const auto ez = square(track_z_error) + square(cluster_z_error);

    // sanity check
    if( std::isnan( erp ) ) continue;
    if( std::isnan( ez ) ) continue;

    // get residuals
    const auto drp = cluster_r*delta_phi( track_phi - cluster_phi );
    const auto dz = track_z - cluster_z;

    // get track angles
    const auto cosphi( std::cos( track_phi ) );
    const auto sinphi( std::sin( track_phi ) );
    const auto track_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto track_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto track_pz = state->get_pz();
    const auto talpha = -track_pphi/track_pr;
    const auto tbeta = -track_pz/track_pr;

    // check against limits
    // TODO: make this configurable
    static constexpr float max_talpha = 0.6;
    static constexpr float max_residual = 5;
    if( std::abs( talpha ) > max_talpha ) continue;
    if( std::abs( drp ) > max_residual ) continue;

    // get cell
    const auto i = get_cell( cluster_key, cluster );

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
void  PHSpaceChargeReconstruction::calculate_distortions( PHCompositeNode* topNode )
{

  // get tpc geometry
  auto *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << " can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return;
  }

  // calculate distortions in each volume elements
  std::vector<column_t> delta(m_totalbins);
  std::vector<matrix_t> cov(m_totalbins);

  for( int i = 0; i < m_totalbins; ++i )
  {
    cov[i] = m_lhs[i].inverse();
    delta[i] = m_lhs[i].partialPivLu().solve( m_rhs[i] );
  }

  // create tgraphs
  using TGraphPointer = std::unique_ptr<TGraphErrors>;
  std::vector<TGraphPointer> tg( m_zbins*m_phibins*m_ncoord );

  for( int iz = 0; iz < m_zbins; ++iz )
    for( int iphi = 0; iphi < m_phibins; ++iphi )
    for( int icoord = 0; icoord < m_ncoord; ++icoord )
  {

    const int tgindex = iz + m_zbins*( iphi + m_phibins*icoord );
    tg[tgindex].reset( new TGraphErrors() );
    tg[tgindex]->SetName( Form( "tg_%i_%i_%i", iz, iphi, icoord ) );

    for( int ir = 0; ir < m_rbins; ++ir )
    {

      // get layers corresponding to bins
      const int inner_layer = m_firstlayer_tpc + m_nlayers_tpc*ir/m_rbins;
      const int outer_layer = m_firstlayer_tpc + m_nlayers_tpc*(ir+1)/m_rbins-1;

      const auto inner_radius = geom_container->GetLayerCellGeom(inner_layer)->get_radius();
      const auto outer_radius = geom_container->GetLayerCellGeom(outer_layer)->get_radius();
      const float r = (inner_radius+outer_radius)/2;

      int index = get_cell( iz, ir, iphi );
      tg[tgindex]->SetPoint( ir, r, delta[index](icoord,0) );
      tg[tgindex]->SetPointError( ir, 0, std::sqrt(cov[index](icoord,icoord)) );
    }
  }

  // save to root file
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    for( auto&& value:tg ) { value->Write(); }
    outputfile->Close();
  }
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::get_cell( int iz, int ir, int iphi ) const
{
  if( ir < 0 || ir >= m_rbins ) return -1;
  if( iphi < 0 || iphi >= m_phibins ) return -1;
  if( iz < 0 || iz >= m_zbins ) return -1;
  return iz + m_zbins*( ir + m_rbins*iphi );
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::get_cell( TrkrDefs::cluskey cluster_key, TrkrCluster* cluster ) const
{
  // radius
  const auto layer = TrkrDefs::getLayer(cluster_key);
  const int ir = m_rbins*(layer - m_firstlayer_tpc)/m_nlayers_tpc;

  // azimuth
  auto cluster_phi = get_phi( cluster->getX(), cluster->getY() );
  if( cluster_phi >= M_PI ) cluster_phi -= 2*M_PI;
  const int iphi = m_phibins*(cluster_phi + M_PI)/(2.*M_PI);

  // z
  const auto cluster_z = cluster->getZ();
  static constexpr float z_min = -212/2;
  static constexpr float z_max = 212/2;
  const int iz = m_zbins*(cluster_z-z_min)/(z_max-z_min);

  return get_cell( iz, ir, iphi );
}

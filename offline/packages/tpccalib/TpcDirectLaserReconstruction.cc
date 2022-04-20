/**
 * \file TpcDirectLaserReconstruction.cc
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDirectLaserReconstruction.h"

#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>

#include <cassert>

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

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! get radius from x and y
  template<class T>
    inline constexpr T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  // calculate intersection between line and circle
  double line_circle_intersection( const TVector3& p, const TVector3& d, double radius )
  {
    const double A = square(d.x()) + square(d.y());
    const double B = 2*p.x()*d.x() + 2*p.y()*d.y();
    const double C = square(p.x()) + square(p.y()) - square(radius);
    const double delta = square(B)-4*A*C;
    if( delta < 0 ) return -1;

    // check first intersection
    const double tup = (-B + std::sqrt(delta))/(2*A);
    if( tup >= 0 ) return tup;

    // check second intersection
    const double tdn = (-B-sqrt(delta))/(2*A);
    if( tdn >= 0 ) return tdn;

    // no valid extrapolation
    return -1;
  }
    
  /// TVector3 stream
  inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

  /// calculate delta_phi between -pi and pi
  template< class T>
    inline constexpr T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
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
TpcDirectLaserReconstruction::TpcDirectLaserReconstruction( const std::string& name ):
  SubsysReco( name)
  , PHParameterInterface(name)
  , m_matrix_container( new TpcSpaceChargeMatrixContainerv1 )
{
  InitializeParameters();
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::Init(PHCompositeNode*)
{
  m_total_clusters = 0;
  m_accepted_clusters = 0;

  if( m_savehistograms ) create_histograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::InitRun(PHCompositeNode* )
{
  UpdateParametersWithMacro();
  m_max_dca = get_double_param( "directlaser_max_dca" );
  m_max_drphi = get_double_param( "directlaser_max_drphi" );
  m_max_dz = get_double_param( "directlaser_max_dz" );

  // print
  if( Verbosity() )
  {
    std::cout
      << "TpcDirectLaserReconstruction::InitRun\n"
      << " m_outputfile: " << m_outputfile << "\n"
      << " m_max_dca: " << m_max_dca << "\n"
      << " m_max_drphi: " << m_max_drphi << "\n"
      << " m_max_dz: " << m_max_dz << "\n"
      << std::endl;

    // also identify the matrix container
    m_matrix_container->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::End(PHCompositeNode* )
{
  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }

  // write evaluation histograms to output
  if( m_savehistograms && m_histogramfile )
  {
    m_histogramfile->cd();
    for(const auto& o:std::initializer_list<TObject*>({ h_dca_layer, h_deltarphi_layer, h_deltaz_layer, h_entries }))
    { if( o ) o->Write(); }
    m_histogramfile->Close();
  }

  // print counters
  std::cout
    << "TpcDirectLaserReconstruction::End -"
    << " cluster statistics total: " << m_total_clusters
    << " accepted: " << m_accepted_clusters << " fraction: "
    << 100.*m_accepted_clusters/m_total_clusters << "%"
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void TpcDirectLaserReconstruction::SetDefaultParameters()
{

  // DCA cut, to decide whether a cluster should be associated to a given laser track or not
  set_default_double_param( "directlaser_max_dca", 1.5 );

  
//   // residual cuts, used to decide if a given cluster is used to fill SC reconstruction matrices
//   set_default_double_param( "directlaser_max_drphi", 0.5 );
//   set_default_double_param( "directlaser_max_dz", 0.5 );

  set_default_double_param( "directlaser_max_drphi", 2. );
  set_default_double_param( "directlaser_max_dz", 2. );
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::set_grid_dimensions( int phibins, int rbins, int zbins )
{ m_matrix_container->set_grid_dimensions( phibins, rbins, zbins ); }

//_____________________________________________________________________
int TpcDirectLaserReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // acts surface map
  m_surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  assert( m_surfmaps );

  // acts geometry
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  assert( m_tGeometry );

  // tracks
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

  // clusters
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert(m_cluster_map);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::create_histograms()
{
  std::cout << "TpcDirectLaserReconstruction::makeHistograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  // residuals vs layers
  h_dca_layer = new TH2F( "dca_layer", ";layer; DCA (cm)", 57, 0, 57, 500, 0, 2 );
  h_deltarphi_layer = new TH2F( "deltarphi_layer", ";layer; r.#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 2000, -2, 2 );
  h_deltaz_layer = new TH2F( "deltaz_layer", ";layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 400, -2, 2 );

  // entries vs cell grid
  /* histogram dimension and axis limits must match that of TpcSpaceChargeMatrixContainer */
  if( m_matrix_container )
  {
    int phibins = 0;
    int rbins = 0;
    int zbins = 0;
    m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );
    h_entries = new TH3F( "entries", ";#phi;r (cm);z (cm)", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax );
  }
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;

  // count number of clusters in the TPC
  for(const auto& hitsetkey:m_cluster_map->getHitSetKeys(TrkrDefs::tpcId))
  {
    const auto range = m_cluster_map->getClusters(hitsetkey);
    m_total_clusters += std::distance( range.first, range.second );
  }

  // loop over tracks and process
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_track( SvtxTrack* track )
{

  // get track parameters
  const TVector3 origin( track->get_x(), track->get_y(), track->get_z() );
  const TVector3 direction( track->get_px(), track->get_py(), track->get_pz() );

  if( Verbosity() )
  { std::cout << "TpcDirectLaserReconstruction::process_track - position: " << origin << " direction: " << direction << std::endl; }

  // loop over hitsets
  for(const auto& hitsetkey:m_cluster_map->getHitSetKeys(TrkrDefs::tpcId))
  {
    
    // store layer
    const auto layer = TrkrDefs::getLayer( hitsetkey );

    // get corresponding clusters
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      
      // get cluster global coordinates
      const auto global = m_transformer.getGlobalPosition(cluster,m_surfmaps, m_tGeometry);

      // calculate dca
      const TVector3 oc( global.x()-origin.x(), global.y()-origin.y(), global.z()-origin.z()  );
      auto t = direction.Dot( oc )/square( direction.Mag() );
      auto om = direction*t;
      const auto dca = (oc-om).Mag();
      
      // do not associate if dca is too large
      if( dca > m_max_dca ) continue;
      
      // calculate intersection to layer
      t = line_circle_intersection(origin, direction, get_r( global.x(), global.y() ));
      if( t < 0 ) continue;

      // update position on track
      om = direction*t;
   
      // path length
      const auto pathlength = om.Mag();

      // get projection to the track
      const auto projection = origin + om;

      // create relevant state vector and assign to track
      SvtxTrackState_v1 state( pathlength );
      state.set_x( projection.x() );
      state.set_y( projection.y() );
      state.set_z( projection.z() );

      state.set_px( direction.x());
      state.set_py( direction.y());
      state.set_pz( direction.z());
      track->insert_state( &state );

      // also associate cluster to track
      track->insert_cluster_key( key );

      // cluster r, phi and z
      const auto cluster_r = get_r(global.x(), global.y());
      const auto cluster_phi = std::atan2(global.y(),global.x());
      const auto cluster_z = global.z();

      // cluster errors
      const auto cluster_rphi_error = cluster->getRPhiError();
      const auto cluster_z_error = cluster->getZError();

//       /*
//       remove clusters with too small errors since they are likely pathological
//       and have a large contribution to the chisquare
//       TODO: make these cuts configurable
//       */
//       if( cluster_rphi_error < 0.015 ) continue;
//       if( cluster_z_error < 0.05 ) continue;

      // track position
      const auto track_phi = std::atan2( projection.y(), projection.x() );
      const auto track_z = projection.z();

      // track angles
      const auto cosphi( std::cos( track_phi ) );
      const auto sinphi( std::sin( track_phi ) );
      const auto track_pphi = -state.get_px()*sinphi + state.get_py()*cosphi;
      const auto track_pr = state.get_px()*cosphi + state.get_py()*sinphi;
      const auto track_pz = state.get_pz();
      const auto talpha = -track_pphi/track_pr;
      const auto tbeta = -track_pz/track_pr;

      // sanity check
      if( std::isnan(talpha) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - talpha is nan" << std::endl;
        continue;
      }

      if( std::isnan(tbeta) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - tbeta is nan" << std::endl;
        continue;
      }

      // residuals
      const auto drp = cluster_r*delta_phi( cluster_phi - track_phi );
      const auto dz = cluster_z - track_z;

      // sanity checks
      if( std::isnan(drp) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - drp is nan" << std::endl;
        continue;
      }

      if( std::isnan(dz) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - dz is nan" << std::endl;
        continue;
      }

      if(m_savehistograms)
      {
        if(h_dca_layer) h_dca_layer->Fill(layer, dca);
        if(h_deltarphi_layer) h_deltarphi_layer->Fill(layer, drp);
        if(h_deltaz_layer) h_deltaz_layer->Fill(layer, dz);
        if(h_entries)
        {
          auto phi = cluster_phi;
          while( phi < m_phimin ) phi += 2.*M_PI;
          while( phi >= m_phimax ) phi -= 2.*M_PI;
          h_entries->Fill( phi, cluster_r, cluster_z );
        }
      }

//       // check against limits
//       if( std::abs( drp ) > m_max_drphi ) continue;
//       if( std::abs( dz ) > m_max_dz ) continue;

      // residual errors squared
      const auto erp = square(cluster_rphi_error);
      const auto ez = square(cluster_z_error);

      // sanity check
      if( std::isnan( erp ) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - erp is nan" << std::endl;
        continue;
      }

      if( std::isnan( ez ) )
      {
        std::cout << "TpcDirectLaserReconstruction::process_track - ez is nan" << std::endl;
        continue;
      }

      // get cell
      const auto i = get_cell_index( global );
      if( i < 0 )
      {
        if( Verbosity() )
        {
          std::cout << "TpcDirectLaserReconstruction::process_track - invalid cell index"
            << " r: " << cluster_r
            << " phi: " << cluster_phi
            << " z: " << cluster_z
            << std::endl;
        }
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

}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::get_cell_index( const Acts::Vector3& global ) const
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );

  // phi
  // bound check
  float phi = std::atan2( global.y(), global.x() );
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  int iphi = phibins*(phi-m_phimin)/(m_phimax-m_phimin);

  // radius
  const float r = get_r( global.x(), global.y() );
  if( r < m_rmin || r >= m_rmax ) return -1;
  int ir = rbins*(r-m_rmin)/(m_rmax-m_rmin);

  // z
  const float z = global.z();
  if( z < m_zmin || z >= m_zmax ) return -1;
  int iz = zbins*(z-m_zmin)/(m_zmax-m_zmin);

  return m_matrix_container->get_cell_index( iphi, ir, iz );
}

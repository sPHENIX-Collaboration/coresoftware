/**
 * \file TpcDirectLaserReconstruction.cc
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDirectLaserReconstruction.h"

#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit
#include <trackbase/TpcDefs.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>
#include <TNtuple.h>

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
  std::pair<double, double> line_circle_intersection( const TVector3& p, const TVector3& d, double radius )
  {
    std::pair<double,double> ret = std::make_pair(-1, -1);

    const double A = square(d.x()) + square(d.y());
    const double B = 2*p.x()*d.x() + 2*p.y()*d.y();
    const double C = square(p.x()) + square(p.y()) - square(radius);
    const double delta = square(B)-4*A*C;
    if( delta < 0 ) return ret;

    // we want 1 intersection
    const double tup = (-B + std::sqrt(delta))/(2*A);
    const double tdn = (-B-sqrt(delta))/(2*A);
    // std::cout << " tup " << tup << " tdn " << tdn << std::endl;
    ret = std::make_pair(tup, tdn);

    return ret;

    /*
    if( tup > 0 && tup < tdn) 
      return tup;
    if(tdn > 0 && tdn < tup)
       return tdn;

    // no valid extrapolation
    return -1;
    */
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
  m_total_hits = 0;
  m_matched_hits = 0;
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
    //for(const auto& o:std::initializer_list<TObject*>({ h_dca_layer, h_deltarphi_layer_south,h_deltarphi_layer_north, h_deltaz_layer, h_deltar_r,h_deltheta_delphi,h_deltheta_delphi_1,h_deltheta_delphi_2,h_deltheta_delphi_3,h_deltheta_delphi_4,h_deltheta_delphi_5,h_deltheta_delphi_6,h_deltheta_delphi_7,h_deltheta_delphi_8, h_entries,h_relangle_lasrangle,h_relangle_theta_lasrangle,h_relangle_phi_lasrangle,h_GEMs_hit,h_layers_hit,h_xy,h_xz,h_xy_pca,h_xz_pca,h_dca_path,h_zr,h_zr_pca, h_dz_z }))

    for(const auto& o:std::initializer_list<TObject*>({ h_dca_layer, h_deltarphi_layer_south,h_deltarphi_layer_north, h_deltaz_layer, h_deltar_r,h_deltheta_delphi,h_deltheta_delphi_1,h_deltheta_delphi_2,h_deltheta_delphi_3,h_deltheta_delphi_4,h_deltheta_delphi_5,h_deltheta_delphi_6,h_deltheta_delphi_7,h_deltheta_delphi_8, h_entries,h_hits,h_hits_reco,h_adc,h_adc_reco,h_adc_sum,h_adc_sum_reco,h_adc_sum_ratio,h_adc_sum_ratio_true,h_adc_vs_DCA_true,h_adc_sum_ratio_lasrangle,h_num_sum,h_num_sum_reco,h_num_sum_ratio,
h_num_sum_ratio_true,h_adc_sum_ratio_true,h_adc_sum_ratio,h_num_sum_ratio_lasrangle,h_bright_hits_laser1,h_bright_hits_laser2,h_bright_hits_laser3,h_bright_hits_laser4,h_relangle_lasrangle,
h_relangle_theta_lasrangle,h_relangle_phi_lasrangle,h_GEMs_hit,h_layers_hit,h_xy,h_xz,h_xy_pca,h_xz_pca,h_dca_path,h_zr,h_zr_pca, h_dz_z }))

    { if( o ) o->Write(); }
    m_histogramfile->Close();
  }

// add these back in later !!!
//,h_assoc_hits
//h_hits
//h_bright_hits_laser1
//h_bright_hits_laser2
//h_bright_hits_laser3
//h_bright_hits_laser4
//h_origins
//h_clusters

  // print counters
  std::cout << "TpcDirectLaserReconstruction::End - m_total_hits: " << m_total_hits << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - m_matched_hits: " << m_matched_hits << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - m_accepted_clusters: " << m_accepted_clusters << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - fraction cluster/hits: " << m_accepted_clusters/m_total_hits << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void TpcDirectLaserReconstruction::SetDefaultParameters()
{

  // DCA cut, to decide whether a cluster should be associated to a given laser track or not (set to 5 - Charles 10/24/23)
  set_default_double_param( "directlaser_max_dca", 1.0 );

  
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
  m_geom_container =  findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert( m_geom_container );

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // tracks
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

 // get node containing the digitized hits
  m_hit_map = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hit_map);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::create_histograms()
{
  std::cout << "TpcDirectLaserReconstruction::makeHistograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  // residuals vs layers
  h_dca_layer = new TH2F( "dca_layer", ";radius; DCA (cm)", 78, 0, 78, 500, 0, 20 );
  h_deltarphi_layer_north = new TH2F( "deltarphi_layer_north", ";radius; r.#Delta#phi_{track-cluster} (cm)", 78, 0, 78, 2000, -5, 5 );
  h_deltarphi_layer_south = new TH2F( "deltarphi_layer_south", ";radius; r.#Delta#phi_{track-cluster} (cm)", 78, 0, 78, 2000, -5, 5 );
  h_deltaz_layer = new TH2F( "deltaz_layer", ";radius; #Deltaz_{track-cluster} (cm)", 78, 0, 78, 2000, -20, 20 );
  h_deltar_r = new TH2F("deltar_r",";radius;#Deltar_{track-cluster} (cm)", 78,0,78,2000,-3,3);

  h_xy = new TH2F("h_xy"," x vs y", 320,-80,80,320,-80,80);
  h_xz = new TH2F("h_xz"," x vs z", 320,-80,80,440,-110,110);
  h_xy_pca = new TH2F("h_xy_pca"," x vs y pca", 320,-80,80,320,-80,80);
  h_xz_pca = new TH2F("h_xz_pca"," x vs z pca", 320,-80,80,440,-110,110);
  h_dca_path = new TH2F("h_dca_path"," dca vs pathlength", 440,0,110,100,0,20);
  h_zr = new TH2F("h_zr"," z vs r", 440,-110,110,1000,28,80);
  h_zr->GetXaxis()->SetTitle("z");
  h_zr->GetYaxis()->SetTitle("rad");
  h_zr_pca = new TH2F("h_zr_pca"," z vs r pca", 440,-110,110,1000,28,80);
  h_dz_z = new TH2F("h_dz_z"," dz vs z", 440,-110,110, 1000, -20, 20);

  h_hits = new TNtuple("hits","raw hits","x:y:z:adc");
  h_hits_reco = new TNtuple("hits_reco","raw hits","x:y:z:adc");

  h_adc = new TH1F("adc","pedestal subtracted ADC spectra of ALL lasers (MCTRUTH direction)",1063,-19.5,1042.5);
  h_adc->SetXTitle("ADC - PEDESTAL [ADU]");

  h_adc_reco = new TH1F("adc_reco","pedestal subtracted ADC spectra of ALL lasers (RECO DIRECTION)",1063,-19.5,1042.5);
  h_adc_reco->SetXTitle("ADC - PEDESTAL [ADU]");

  h_adc_sum = new TH1F("adc_sum","non-pedestal subtracted sum of ADC for each laser (MCTRUTH DIRECTION)",1600,-0.5,159999.5);
  h_adc_sum->SetXTitle("#Sigma ADC [ADU]");

  h_adc_sum_reco = new TH1F("adc_sum_reco","non-pedestal subtracted sum of ADC for each laser (RECO DIRECTION)",1600,-0.5,159999.5);
  h_adc_sum_reco->SetXTitle("#Sigma ADC [ADU]");

  h_num_sum = new TH1F("num_sum","sum of hits for each laser (MCTRUTH DIRECTION)",1600,-0.5,7999.5);
  h_num_sum->SetXTitle("#Sigma N_{hits}^{laser X} [arb. units]");

  h_num_sum_reco = new TH1F("num_sum_reco","sum of hits for each laser (RECO DIRECTION)",1600,-0.5,7999.5);
  h_num_sum_reco->SetXTitle("#Sigma N_{hits}^{laser X} [arb. units]");

  h_adc_sum_ratio = new TH1F("adc_sum_ratio","RATIO of non-pedestal subtracted sum of ADC for each laser (RECO DIRECTION/MCTRUTH DIRECTION)",120,-0.5,5.5);
  h_adc_sum_ratio->SetXTitle("#frac{#Sigma^{RECO} ADC}{#Sigma^{MCTRUTH} ADC}  [arb. units]");

  h_adc_sum_ratio_true = new TH1F("adc_sum_ratio_true","RATIO of non-pedestal subtracted sum of ADC for each laser (MCTRUTH DIRECTION/ALL )",120,-0.5,5.5);
  h_adc_sum_ratio_true->SetXTitle("#frac{#Sigma^{MCTRUTH} ADC}{#Sigma^{ALL} ADC}  [arb. units]");

  h_num_sum_ratio = new TH1F("num_sum_ratio","RATIO of number of hits associated to each laser (RECO DIRECTION/MCTRUTH DIRECTION)",120,-0.5,5.5);
  h_num_sum_ratio->SetXTitle("#frac{#Sigma^{RECO} N_{hits}^{laser X}}{#Sigma^{MCTRUTH} N_{hits}^{laser X}}  [arb. units]");

  h_num_sum_ratio_true = new TH1F("num_sum_ratio_true","RATIO of number of hits associated to each laser (MCTRUTH DIRECTION/ALL )",120,-0.5,5.5);
  h_num_sum_ratio_true->SetXTitle("#frac{#Sigma^{MCTRUTH} N_{hits}^{laser X}}{#Sigma^{ALL} N_{hits}^{laser X}}  [arb. units]");

  h_adc_vs_DCA_true = new TH2F("adc_vs_DCA_true","PEDESTAL SUBTRACTED ADC vs DCA for ALL HITS from ALL LASERS",100,0,100,1024,-0.5,1023.5);
  h_adc_vs_DCA_true->SetXTitle("DCA [mm]");
  h_adc_vs_DCA_true->SetYTitle("ADC - PEDESTAL [ADU]");

  h_adc_sum_ratio_lasrangle = new TH2F("adc_sum_ratio_lasrangle","RATIO of non-pedestal subtracted sum of ADC for LASER 1 ONLY (RECO DIRECTION/MCTRUTH DIRECTION)",91,-0.5,90.5,361,-0.5,360.5);
  h_adc_sum_ratio_lasrangle->SetXTitle("#theta_{laser} (deg.)");
  h_adc_sum_ratio_lasrangle->SetYTitle("#phi_{laser} (deg.)");

  h_num_sum_ratio_lasrangle = new TH2F("num_sum_ratio_lasrangle","RATIO of sum of hits for LASER 1 ONLY (RECO DIRECTION/MCTRUTH DIRECTION)",91,-0.5,90.5,361,-0.5,360.5);
  h_num_sum_ratio_lasrangle->SetXTitle("#theta_{laser} (deg.)");
  h_num_sum_ratio_lasrangle->SetYTitle("#phi_{laser} (deg.)");

  h_bright_hits_laser1 = new TNtuple("bright_hits_laser1","bright hits relative to laser 1","x:y:z:deltheta:delphi");
  h_bright_hits_laser2 = new TNtuple("bright_hits_laser2","bright hits relative to laser 2","x:y:z:deltheta:delphi");
  h_bright_hits_laser3 = new TNtuple("bright_hits_laser3","bright hits relative to laser 3","x:y:z:deltheta:delphi");
  h_bright_hits_laser4 = new TNtuple("bright_hits_laser4","bright hits relative to laser 4","x:y:z:deltheta:delphi");

  //h_assoc_hits = new TNtuple("assoc_hits","hits associated with tracks (dca cut)","x:y:z");
  //h_clusters = new TNtuple("clusters","associated clusters","x:y:z");
  //h_origins = new TNtuple("origins","track origins","x:y:z");

  //h_relangle_lasrangle = new TH2F("relangle_lasrangle","Difference in pointing angle from laser and relative angle histogram for ALL lasers",51,-25.5,25.5,51,-25.5,25.5);
  h_relangle_lasrangle = new TH2F("relangle_lasrangle","Difference in pointing angle from laser and relative angle histogram for ALL lasers",102,-25.5,25.5,102,-25.5,25.5);
  h_relangle_lasrangle->SetXTitle("#delta#theta (deg.)");
  h_relangle_lasrangle->SetYTitle("#delta#phi (deg.)");

  h_relangle_theta_lasrangle = new TH2F("relangle_theta_lasrangle","#delta#theta from laser and relative angle histogram for LASER 1 ONLY, ",91,-0.5,90.5,361,-0.5,360.5);
  h_relangle_theta_lasrangle->SetXTitle("#theta_{laser} (deg.)");
  h_relangle_theta_lasrangle->SetYTitle("#phi_{laser} (deg.)");

  h_relangle_phi_lasrangle = new TH2F("relangle_phi_lasrangle","#delta#phi from laser and relative angle histogram for LASER 1 ONLY",91,-0.5,90.5,361,-0.5,360.5);
  h_relangle_phi_lasrangle->SetXTitle("#theta_{laser} (deg.)");
  h_relangle_phi_lasrangle->SetYTitle("#phi_{laser} (deg.)");

  h_deltheta_delphi = new TH2F("deltheta_delphi","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and ALL laser start points", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi->SetXTitle("#Delta#theta");
  h_deltheta_delphi->SetYTitle("#Delta#phi");

  h_deltheta_delphi_1 = new TH2F("deltheta_delphi_1","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 1 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_1->SetXTitle("#Delta#theta");
  h_deltheta_delphi_1->SetYTitle("#Delta#phi");

  h_deltheta_delphi_2 = new TH2F("deltheta_delphi_2","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 2 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_2->SetXTitle("#Delta#theta");
  h_deltheta_delphi_2->SetYTitle("#Delta#phi");

  h_deltheta_delphi_3 = new TH2F("deltheta_delphi_3","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 3 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_3->SetXTitle("#Delta#theta");
  h_deltheta_delphi_3->SetYTitle("#Delta#phi");

  h_deltheta_delphi_4 = new TH2F("deltheta_delphi_4","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 4 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_4->SetXTitle("#Delta#theta");
  h_deltheta_delphi_4->SetYTitle("#Delta#phi");

  h_deltheta_delphi_5 = new TH2F("deltheta_delphi_5","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 5 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_5->SetXTitle("#Delta#theta");
  h_deltheta_delphi_5->SetYTitle("#Delta#phi");

  h_deltheta_delphi_6 = new TH2F("deltheta_delphi_6","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 6 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_6->SetXTitle("#Delta#theta");
  h_deltheta_delphi_6->SetYTitle("#Delta#phi");

  h_deltheta_delphi_7 = new TH2F("deltheta_delphi_7","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 7 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_7->SetXTitle("#Delta#theta");
  h_deltheta_delphi_7->SetYTitle("#Delta#phi");

  h_deltheta_delphi_8 = new TH2F("deltheta_delphi_8","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 8 only", 181, -0.5,180.5,361,-0.5,360.5);
  h_deltheta_delphi_8->SetXTitle("#Delta#theta");
  h_deltheta_delphi_8->SetYTitle("#Delta#phi");


  h_GEMs_hit = new TH1F("GEMS_hit","Number of Unique GEM Modules hit for each laser",8,0,8);
  h_layers_hit = new TH1F("layers_hit","Number of Unique Layers hit for each laser",8,8,8);

  char GEM_bin_label[128];
  for(int GEMhistiter = 0; GEMhistiter < 8; GEMhistiter++)
  {  // (pos z) laser 1 {0,60}, laser 2 {60,0}, laser 3 {0,-60}, laser 4 {-60,0}, (neg z) laser 5 {0,60}, laser 2 {60,0}, laser 3 {0,-60}, laser 4 {-60,0}
    sprintf(GEM_bin_label,"laser %i",GEMhistiter+1);
    h_GEMs_hit->GetXaxis()->SetBinLabel(GEMhistiter+1,GEM_bin_label);
    h_layers_hit->GetXaxis()->SetBinLabel(GEMhistiter+1,GEM_bin_label);
  }
  h_GEMs_hit->SetYTitle("Number of Unique GEM Modules Hit");
  h_layers_hit->SetYTitle("Number of Unique Layers Hit");

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
  if( !( m_track_map && m_hit_map ) ) return;

  // loop over tracks and process
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_track( SvtxTrack* track )
{
  std::multimap<unsigned int, std::pair<float, TVector3>> cluspos_map;
  std::set<unsigned int> layer_bin_set;

  // get track parameters
  const TVector3 origin( track->get_x(), track->get_y(), track->get_z() );
  const TVector3 direction( track->get_px(), track->get_py(), track->get_pz() );
  TVector3 direction2( track->get_px(), track->get_py(), track->get_pz() );

/*
  if(h_origins)
  {
    h_origins->Fill(track->get_x(), track->get_y(), track->get_z());
  }
*/
  const unsigned int trkid = track->get_id();

  const float theta_orig = origin.Theta();
  const float phi_orig = origin.Phi();

  float theta_trk = direction.Theta();
  direction2.RotateZ(-phi_orig);
  float   phi_trk = direction2.Phi();

  if(theta_trk > M_PI/2.) theta_trk = M_PI - theta_trk;

  //if(  phi_trk < 0) phi_trk*=-1;

  while( phi_trk < m_phimin ) phi_trk += 2.*M_PI; 
  while( phi_trk >= m_phimax ) phi_trk -= 2.*M_PI; 

  //const double xo = track->get_x();
  //const double yo = track->get_y();
  //const double zo = track->get_z();

  //if( phi_orig < 0 )phi_orig += 2.*M_PI;

  if( track->get_id() > 3 ){ phi_trk = (2 * M_PI) - phi_trk;}

  //if( Verbosity() )
  //{ 
    std::cout << "----------  processing track " << track->get_id() << std::endl;
    std::cout << "TpcDirectLaserReconstruction::process_track - position: " << origin << " direction: " << direction << std::endl; 
    std::cout << "Position Orientation - Theta: " << theta_orig*(180./M_PI) << ", Phi: " << phi_orig*(180./M_PI) << std::endl;
    std::cout << "Track Angle Direction - Theta: " << theta_trk*(180./M_PI) << ", Phi: " << phi_trk*(180./M_PI) << std::endl;
  //}
  
  int GEM_Mod_Arr[72] = {0};

  // loop over hits
  TrkrHitSetContainer::ConstRange hitsetrange = m_hit_map->getHitSets(TrkrDefs::TrkrId::tpcId);

  float max_adc = 0.;
  float sum_adc_truth=0;
  float sum_adc_truth_all=0;

  float sum_n_hits_truth=0;
  float sum_n_hits_truth_all=0;

  for(auto hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    const TrkrDefs::hitsetkey& hitsetkey = hitsetitr->first;
    const int side = TpcDefs::getSide(hitsetkey);
    
    auto hitset = hitsetitr->second;
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    const auto layergeom = m_geom_container->GetLayerCellGeom(layer);
    const auto layer_center_radius = layergeom->get_radius();
    
    // maximum drift time.
    /* it is needed to calculate a given hit position from its drift time */
    static constexpr double AdcClockPeriod = 53.0;   // ns 
    const unsigned short NTBins = (unsigned short)layergeom->get_zbins();
    const float tdriftmax =  AdcClockPeriod * NTBins / 2.0;

    // get corresponding hits
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    
    for (auto hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
      {
	++m_total_hits;
        sum_n_hits_truth_all++;

	const unsigned short phibin = TpcDefs::getPad(hitr->first);
  const unsigned short zbin = TpcDefs::getTBin(hitr->first);

	const double phi = layergeom->get_phicenter(phibin);
  const double x = layer_center_radius * cos(phi);
  const double y = layer_center_radius * sin(phi);
  
  const double zdriftlength = layergeom->get_zcenter(zbin)*m_tGeometry->get_drift_velocity();
  double z  =  tdriftmax*m_tGeometry->get_drift_velocity() - zdriftlength;
  if(side == 0)  z *= -1;
      
	const TVector3 global(x,y,z);

	float adc = (hitr->second->getAdc()) - m_pedestal;
        float adc_unsub = (hitr->second->getAdc()); 
        sum_adc_truth_all += adc_unsub;
/*
        if(h_hits && trkid == 0)
          {
            h_hits->Fill(x,y,z,adc);
          }
*/
        if( adc > max_adc ){ max_adc = adc; }
	
	// calculate dca
	// origin is track origin, direction is track direction
        bool sameside = sameSign(global.z(),origin.z());
	const TVector3 oc( global.x()-origin.x(), global.y()-origin.y(), global.z()-origin.z()  );  // vector from track origin to cluster
	auto t = direction.Dot( oc )/square( direction.Mag() );
	auto om = direction*t;     // vector from track origin to PCA
	const auto dca = (oc-om).Mag();

        if (sameside){ h_adc_vs_DCA_true->Fill(dca,adc);}

        //rotate displacement vector by the rotation of the coordinates:
        //oc2.RotateY(-theta_orig * trkzdir); //undoing rotation in PHG4TpcDirectLaser 

        //global2.RotateZ(-phi_orig);
        //origin2.RotateZ(-phi_orig);

        TVector3 oc2( global.x()-origin.x(), global.y()-origin.y(), global.z()-origin.z()  );  // vector from track origin to cluster
        //const double xt = x-xo;
        //const double yt = y-yo;
        //const double zt = z-zo;
        oc2.RotateZ(-phi_orig);
   

        //auto theta1 = origin.Theta()*(180./M_PI);
        //auto theta2 = global.Theta()*(180./M_PI);

        //auto phi1 = origin.Phi()*(180./M_PI);
        //auto phi2 = global.Phi()*(180./M_PI);

        //float deltheta = GetRelTheta( origin.x() , origin.y(), origin.z(), global.x() , global.y() , global.z() , origin.Theta(), origin.Phi()) *(180./M_PI); 
        //float delphi = GetRelPhi(origin.x() , origin.y() , global.x() , global.y() , origin.Phi()) *(180./M_PI);

        float deltheta = oc2.Theta();
        float delphi = oc2.Phi();

        if (deltheta > M_PI/2.) deltheta = M_PI - deltheta;

        while( delphi < m_phimin ) delphi += 2.*M_PI; 
        while( delphi >= m_phimax ) delphi -= 2.*M_PI; 

        if( track->get_id() > 3 ){ delphi = (2 * M_PI) - delphi;}

        //if (delphi < 0) delphi*=-1;  

        deltheta *= (180./M_PI);
        delphi   *= (180./M_PI);

/*
        if( trkid < 4)
        {
          deltheta = 180 - theta2;    
          if ( trkid == 1) delphi = fabs(phi2 - phi);
          else delphi = phi2 - phi1;
        }
        else
        {
          deltheta = theta2;
          if (trkid == 7 ) delphi = 360. - fabs(phi2 - phi1);
          else delphi = fabs(phi2 - phi1);
        }
*/

        // relative angle histogram - only fill for hits in the same side as origin

        if(sameside){
         
          //std::cout<<"Trk ID = "<<trkid<<std::endl;

          h_deltheta_delphi->Fill(deltheta, delphi);

          if(track->get_id() == 0) //only fill for 1 laser !!!
          {
            //std::cout<<"Trk ID = "<<trkid<<std::endl; 
            h_deltheta_delphi_1->Fill(deltheta, delphi) ;
/*
            //if( h_bright_hits && deltheta > 80 && deltheta < 95 && delphi > 130 && delphi < 150)
            //filling bright_hits for theta = 75, phi = 80 (simple debugging only - to be removed)
            if( h_bright_hits_laser1 )
            {   
              h_bright_hits_laser1->Fill(x,y,z,deltheta,delphi);
            } 
*/
          }
          if(trkid == 1) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_2->Fill(deltheta, delphi) ;
/*
            if( h_bright_hits_laser2 )
            {   
              h_bright_hits_laser2->Fill(x,y,z,deltheta,delphi);
            } 
*/
          }
          if(trkid == 2) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_3->Fill(deltheta, delphi) ;
/*
            if( h_bright_hits_laser3 )
            {   
              h_bright_hits_laser3->Fill(x,y,z,deltheta,delphi);
            } 
*/
          }
          if(trkid == 3) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_4->Fill(deltheta, delphi) ;
/*
            if( h_bright_hits_laser4 )
            {   
              h_bright_hits_laser4->Fill(x,y,z,deltheta,delphi);
            } 
*/
          }
          if(trkid == 4) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_5->Fill(deltheta, delphi) ;
          }
          if(trkid == 5) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_6->Fill(deltheta, delphi) ;
          }
          if(trkid == 6) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_7->Fill(deltheta, delphi) ;
          }
          if(trkid == 7) //only fill for 1 laser !!!
          {
            h_deltheta_delphi_8->Fill(deltheta, delphi) ;
          }

        }

	// do not associate if dca is too large
	if( dca > m_max_dca ) continue; 

        if(sameside){sum_n_hits_truth++;h_adc->Fill(adc);sum_adc_truth += adc_unsub;} //increment truth BUT ONLY for hits on the same side as origin CHARLES 10.31.23

        ++m_matched_hits;
/*
        if(h_assoc_hits){
          h_assoc_hits->Fill( x,y,z);
        }
*/
        //for locating the associated hits   
        float r2 = std::sqrt((x*x) + (y*y));
        float phi3 = phi;
        float z2 = z;
        while( phi3 < m_phimin ) phi3 += 2.*M_PI; 
        while( phi3 >= m_phimax ) phi3 -= 2.*M_PI; 

        const int locateid = Locate(r2 , phi3 , z2); //find where the cluster is 
      
        if( z2 > m_zmin || z2 < m_zmax ) GEM_Mod_Arr[locateid-1]++; // the array ath the cluster location - counts the number of clusters in each array, (IFF its in the volume !!!)

	// bin hits by layer
	const auto cluspos_pair = std::make_pair(adc, global); 	
	cluspos_map.insert(std::make_pair(layer, cluspos_pair));	
	layer_bin_set.insert(layer);
      } //end looping over hits
  } //end looping over hitset

  h_adc_sum->Fill(sum_adc_truth);
  h_num_sum->Fill(sum_n_hits_truth);

  if(sum_adc_truth_all > 0) // you MUST have hits in this to take ratio
  {
    h_adc_sum_ratio_true->Fill(sum_adc_truth/sum_adc_truth_all); //for a hits associated with the same laser fire take the ratio
  }

  if(sum_n_hits_truth_all > 0) // you MUST have hits in this to take ratio
  {
    h_num_sum_ratio_true->Fill(sum_n_hits_truth/sum_n_hits_truth_all); //for a hits associated with the same laser fire take the ratio
  }

  int maxbin;
  int deltheta_max, delphi_max, dummy_z;

  float theta_reco =0;
  float phi_reco =0;

  float direction_reco;
  if( trkid < 4 ){direction_reco = -1;} //first four lasers are on positive z readout plane, and shoot towards negative z
  else{direction_reco = 1;} // next four lasers are on negative z readout plane and shoot towards positive z

  TVector3 dir(0, 0, direction_reco); //unit vector starts out pointing along z axis

  if( (trkid == 0 && h_deltheta_delphi_1->GetEntries() > 0) && max_adc > 10)
  {
    h_deltheta_delphi_1->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_1->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_1->GetMaximumBin();

    h_deltheta_delphi_1->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" naieve deltheta max: "<< h_deltheta_delphi_1->GetXaxis()->GetBinCenter(deltheta_max) <<", naieve delphi max: "<< h_deltheta_delphi_1->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_1->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_1->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_1->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_1->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_1->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_1->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 1 && h_deltheta_delphi_2->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_2->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_2->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_2->GetMaximumBin();

    h_deltheta_delphi_2->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" naieve deltheta max: "<< h_deltheta_delphi_2->GetXaxis()->GetBinCenter(deltheta_max) <<", naive delphi max: "<< h_deltheta_delphi_2->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_2->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_2->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_2->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_2->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_2->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_2->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 2 && h_deltheta_delphi_3->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_3->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_3->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_3->GetMaximumBin();

    h_deltheta_delphi_3->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" naieve deltheta max: "<< h_deltheta_delphi_3->GetXaxis()->GetBinCenter(deltheta_max) <<", naive delphi max: "<< h_deltheta_delphi_3->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_3->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_3->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_3->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_3->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_3->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_3->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 3 && h_deltheta_delphi_4->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_4->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_4->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_4->GetMaximumBin();

    h_deltheta_delphi_4->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" naieve deltheta max: "<< h_deltheta_delphi_4->GetXaxis()->GetBinCenter(deltheta_max) <<", naieve delphi max: "<< h_deltheta_delphi_4->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_4->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_4->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_4->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_4->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_4->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_4->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 4 && h_deltheta_delphi_5->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_5->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_5->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_5->GetMaximumBin();

    h_deltheta_delphi_5->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" deltheta max: "<< h_deltheta_delphi_5->GetXaxis()->GetBinCenter(deltheta_max) <<", delphi max: "<< h_deltheta_delphi_5->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_5->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_5->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_5->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_5->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_5->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_5->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 5 && h_deltheta_delphi_6->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_6->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_6->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_6->GetMaximumBin();

    h_deltheta_delphi_6->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" deltheta max: "<< h_deltheta_delphi_6->GetXaxis()->GetBinCenter(deltheta_max) <<", delphi max: "<< h_deltheta_delphi_6->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_6->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_6->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_6->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_6->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_6->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_6->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 6 && h_deltheta_delphi_7->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_7->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_7->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_7->GetMaximumBin();

    h_deltheta_delphi_7->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" deltheta max: "<< h_deltheta_delphi_7->GetXaxis()->GetBinCenter(deltheta_max) <<", delphi max: "<< h_deltheta_delphi_7->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_7->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_7->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_7->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_7->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_7->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_7->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  else if( (trkid == 7 && h_deltheta_delphi_8->GetEntries() > 0) && max_adc > 10 )
  {
    h_deltheta_delphi_8->GetXaxis()->SetRangeUser(-5+(theta_trk*(180./M_PI)),5+(theta_trk*(180./M_PI)));
    h_deltheta_delphi_8->GetYaxis()->SetRangeUser(-5+(phi_trk*(180./M_PI)),5+(phi_trk*(180./M_PI)));

    maxbin= h_deltheta_delphi_8->GetMaximumBin();

    h_deltheta_delphi_8->GetBinXYZ(maxbin, deltheta_max, delphi_max, dummy_z);

    std::cout << "Laser # "<<(trkid+1)<<" deltheta max: "<< h_deltheta_delphi_8->GetXaxis()->GetBinCenter(deltheta_max) <<", delphi max: "<< h_deltheta_delphi_8->GetYaxis()->GetBinCenter(delphi_max) << std::endl;

    h_relangle_lasrangle->Fill(h_deltheta_delphi_8->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)),h_deltheta_delphi_8->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)));

    //h_relangle_theta_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_8->GetXaxis()->GetBinCenter(deltheta_max)-(theta_trk*(180./M_PI)) );
    //h_relangle_phi_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),h_deltheta_delphi_8->GetYaxis()->GetBinCenter(delphi_max)-(phi_trk*(180./M_PI)) );

    theta_reco = h_deltheta_delphi_8->GetXaxis()->GetBinCenter(deltheta_max) * (M_PI/180); //convert to radians
    phi_reco = h_deltheta_delphi_8->GetYaxis()->GetBinCenter(delphi_max) * (M_PI/180); //convert to radians
  }

  // Determining Vector for reconstructed laser direction
  dir.RotateY(theta_reco* direction_reco);
  if(direction_reco == -1) dir.RotateZ(phi_reco); //if +z facing -z
  else dir.RotateZ(-phi_reco); //if -z facting +z
  dir.RotateZ(phi_orig); //rotate by the laser coordinate

  //////    print out the reconstructed track direction ////////////////////////////////////////////

  TVector3 direction2_reco( dir.Px(), dir.Py(), dir.Pz() );

  float theta_trk_reco = dir.Theta();
  direction2_reco.RotateZ(-phi_orig);
  float   phi_trk_reco = direction2_reco.Phi();

  if(theta_trk_reco > M_PI/2.) theta_trk_reco = M_PI - theta_trk_reco;

  //if(  phi_trk < 0) phi_trk*=-1;

  while( phi_trk_reco < m_phimin ) phi_trk_reco += 2.*M_PI; 
  while( phi_trk_reco >= m_phimax ) phi_trk_reco -= 2.*M_PI; 

  //const double xo = track->get_x();
  //const double yo = track->get_y();
  //const double zo = track->get_z();

  //if( phi_orig < 0 )phi_orig += 2.*M_PI;

  if( track->get_id() > 3 ){ phi_trk_reco = (2 * M_PI) - phi_trk_reco;}

  std::cout << "RECONSTRUCTED Track Angle Direction - Theta: " << theta_trk_reco*(180./M_PI) << ", Phi: " << phi_trk_reco*(180./M_PI) << std::endl;

  ////////////////////////////////////////////////

  //now loop over hits AGAIN and get ADC spectrum

  float sum_adc_reco=0;
  float sum_n_hits_reco=0;

  for(auto hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    const TrkrDefs::hitsetkey& hitsetkey_2 = hitsetitr->first;
    const int side_2 = TpcDefs::getSide(hitsetkey_2);
    
    auto hitset_2 = hitsetitr->second;
    const unsigned int layer_2 = TrkrDefs::getLayer(hitsetkey_2);
    const auto layergeom_2 = m_geom_container->GetLayerCellGeom(layer_2);
    const auto layer_center_radius_2 = layergeom_2->get_radius();
    
    // maximum drift time.
    /* it is needed to calculate a given hit position from its drift time */
    static constexpr double AdcClockPeriod_2 = 53.0;   // ns 
    const unsigned short NTBins_2 = (unsigned short)layergeom_2->get_zbins();
    const float tdriftmax_2 =  AdcClockPeriod_2 * NTBins_2 / 2.0;

    // get corresponding hits
    TrkrHitSet::ConstRange hitrangei_2 = hitset_2->getHits();
    
    for (auto hitr = hitrangei_2.first; hitr != hitrangei_2.second; ++hitr)
    {

	const unsigned short phibin_2 = TpcDefs::getPad(hitr->first);
        const unsigned short zbin_2 = TpcDefs::getTBin(hitr->first);

	const double phi_2 = layergeom_2->get_phicenter(phibin_2);
        const double x_2 = layer_center_radius_2 * cos(phi_2);
        const double y_2 = layer_center_radius_2 * sin(phi_2);
  
        const double zdriftlength_2 = layergeom_2->get_zcenter(zbin_2)*m_tGeometry->get_drift_velocity();
        double z_2  =  tdriftmax_2*m_tGeometry->get_drift_velocity() - zdriftlength_2;
        if(side_2 == 0)  z_2 *= -1;
      
	const TVector3 global_2(x_2,y_2,z_2);

        bool sameside_reco = sameSign(global_2.z(),origin.z());

	float adc_2 = (hitr->second->getAdc()) - m_pedestal;
        float adc_2_unsub = (hitr->second->getAdc());
	
	// calculate dca
	// origin is track origin, direction is track direction
	const TVector3 oc_2( global_2.x()-origin.x(), global_2.y()-origin.y(), global_2.z()-origin.z()  );  // vector from track origin to cluster
	auto t_2 = dir.Dot( oc_2 )/square( dir.Mag() );
	auto om_2 = dir*t_2;     // vector from track origin to PCA
	const auto dca_2 = (oc_2-om_2).Mag();

        if( dca_2 > m_max_dca ) continue; // do not record if hit is outside DCA, proceed to next hit
        if(sameside_reco){sum_n_hits_reco++;h_adc_reco->Fill(adc_2);sum_adc_reco += adc_2_unsub;} //increment reco but only if hits are on the same side
/*
        if(h_hits_reco && trkid == 0)
        {
            h_hits_reco->Fill(x_2,y_2,z_2,adc_2);
        }
*/
    } //end loop over hits again
  } //end loop over hitset again

  h_adc_sum_reco->Fill(sum_adc_reco);
  h_num_sum_reco->Fill(sum_n_hits_reco);

  if(sum_adc_truth > 0) // you MUST have hits in this to take ratio
  {
    h_adc_sum_ratio->Fill(sum_adc_reco/sum_adc_truth); //for a hits associated with the same laser fire take the ratio
    if( trkid == 0 ) //for laser 1 only, show me WHERE we are messing up reconstruction
    {
      h_adc_sum_ratio_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),(sum_adc_reco/sum_adc_truth));
    }
  }

  if(sum_n_hits_truth > 0) // you MUST have hits in this to take ratio
  {
    h_num_sum_ratio->Fill(sum_n_hits_reco/sum_n_hits_truth); //for a hits associated with the same laser fire take the ratio
    if( trkid == 0 ) //for laser 1 only, show me WHERE we are messing up reconstruction
    {
      h_num_sum_ratio_lasrangle->Fill(theta_trk*(180./M_PI),phi_trk*(180./M_PI),(sum_n_hits_reco/sum_n_hits_truth));
    }
  }

  for(int GEMS_iter = 0; GEMS_iter < 72; GEMS_iter++)
  {
    if(GEM_Mod_Arr[GEMS_iter] > 0){  //laser 0
      h_GEMs_hit->Fill(trkid + 0.5);
    }
  }

  h_deltheta_delphi_1->Reset("ICESM");
  h_deltheta_delphi_2->Reset("ICESM");
  h_deltheta_delphi_3->Reset("ICESM");
  h_deltheta_delphi_4->Reset("ICESM");
  h_deltheta_delphi_5->Reset("ICESM");
  h_deltheta_delphi_6->Reset("ICESM");
  h_deltheta_delphi_7->Reset("ICESM");
  h_deltheta_delphi_8->Reset("ICESM");


  
  // We will bring this back here when we fix the clustering
  //int GEM_Mod_Arr[72] = {0};

  // we have all of the residuals for this track added to the map, binned by layer
  // now we calculate the centroid of the clusters in each bin

  for(auto layer : layer_bin_set)
    {
      PHG4TpcCylinderGeom *layergeom = m_geom_container->GetLayerCellGeom(layer);
      const auto layer_center_radius = layergeom->get_radius();
      const auto layer_inner_radius = layer_center_radius - layergeom->get_thickness() / 2.0;
      const auto layer_outer_radius = layer_center_radius + layergeom->get_thickness() / 2.0;

      h_layers_hit->Fill(trkid + 0.5);
          
      // does track pass completely through this layer? If not do not use the hits
      auto tupdn = line_circle_intersection(origin, direction, layer_outer_radius);
      if( tupdn.first <= 0 && tupdn.second <= 0)
	{
	  std::cout << " punt:  layer " << layer << " layer outer radius " << layer_outer_radius << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}
      double layer_entry = tupdn.first;
      if(tupdn.second >= 0 && tupdn.second < tupdn.first) layer_entry = tupdn.second; 

      tupdn = line_circle_intersection(origin, direction, layer_inner_radius);
      if( tupdn.first <= 0 && tupdn.second <= 0)
	{
	  std::cout << " punt:  layer " << layer << " layer inner radius " << layer_inner_radius << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}
      double layer_exit = tupdn.first;
      if(tupdn.second > 0 && tupdn.second < tupdn.first) layer_exit = tupdn.second;
      if(Verbosity() > 2) std::cout << " layer " << layer << " layer entry " << layer_entry << " layer exit " << layer_exit << std::endl;     

      // calculate track intersection with layer center
      tupdn = line_circle_intersection(origin, direction, layer_center_radius);
      double tup = tupdn.first;
      double tdn = tupdn.second;
      // take 1 one if there are two intersections
      double t = tup;
      if(tdn > 0 && tdn < tup) t = tdn;
      if( t < 0 )
	{
	  std::cout << " punt:  layer " << layer << " layer center radius " << layer_center_radius << " t " << t << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}

      // get position of layer center on track
      auto om = direction*t;      

      // get vector pointing to the track intersection with layer center
      const auto projection = origin + om;
      
      double zmax = -999.;
      double zmin = 999.;
      
      auto bin_residual = cluspos_map.equal_range(layer);

      TVector3 clus_centroid(0,0,0);
      float wt = 0;
      for( auto mit = bin_residual.first; mit != bin_residual.second; ++mit)
	{
	  float adc =  (float) mit->second.first;
	  TVector3 cluspos =mit->second.second;

	  if(fabs(cluspos.z() - projection.z()) > m_max_zrange) continue;  // to reject clusters from possible second traverse of layer
 
	  if(Verbosity() > 2)
	    {
	      std::cout << "  layer " << layer << " adc " << adc << std::endl;
	      std::cout << "            cluspos " << cluspos.x() << "  " << cluspos.y() << "  " << cluspos.z() 
			<< " clus radius " << sqrt(square(cluspos.x()) + square(cluspos.y())) << std::endl;
	    }

	  clus_centroid += cluspos * adc;
	  wt += adc;

	  if(cluspos.z() < zmin) zmin = cluspos.z();
	  if(cluspos.z() > zmax) zmax = cluspos.z();
	}

      clus_centroid.SetX(clus_centroid.x()/ wt);
      clus_centroid.SetY(clus_centroid.y()/ wt);
      clus_centroid.SetZ(clus_centroid.z()/ wt);

      double zrange =  zmax-zmin;
      if(zrange > m_max_zrange)
	{
	  std::cout << "    exceeded  max zrange:  zrange " << zrange << " max zrange " << m_max_zrange << std::endl; 
	  continue;
	}

      // get the distance of the cluster centroid to the track-layer intersection point

      const TVector3 oc( clus_centroid.x()-origin.x(), clus_centroid.y()-origin.y(), clus_centroid.z()-origin.z()  );  // vector from track origin to cluster

      const auto dca = (oc-om).Mag();

      // path length
      const auto pathlength = om.Mag();

      // Correct cluster z for the track transit time using the pathlength 
      double ns_per_cm = 1e9 / 3e10;
      double dt = pathlength * ns_per_cm;
      double transit_dz = dt * m_tGeometry->get_drift_velocity();
      if(origin.z() > 0)
	clus_centroid.SetZ(clus_centroid.z() + transit_dz);
      else
	clus_centroid.SetZ(clus_centroid.z() - transit_dz);

      if(Verbosity() > 0)
	{
	  std::cout << "  layer " << layer << " radius " << layer_center_radius << " wt " << wt << " t " << t << " dca " << dca 
		    << " pathlength " << pathlength << " transit_dz " << transit_dz << std::endl;
	  std::cout << "      clus_centroid " << clus_centroid.x() << "  " << clus_centroid.y() << "  " << clus_centroid.z() << " zrange " << zrange << std::endl;
	  std::cout << "      projection " << projection.x() << "  " << projection.y() << "  " << projection.z() << " dz " << clus_centroid.z() - projection.z() << std::endl;
	}

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
      //track->insert_cluster_key( key );
      
      // cluster r, phi and z
      const auto cluster_r = get_r(clus_centroid.x(), clus_centroid.y());
      const auto cluster_phi = std::atan2(clus_centroid.y(),clus_centroid.x());
      const auto cluster_z = clus_centroid.z();

      // We will bring this back when we fix clustering !!
      /*
      float r2 = cluster_r;
      float phi2 = cluster_phi;
      float z2 = cluster_z;
      while( phi2 < m_phimin ) phi2 += 2.*M_PI; 
      while( phi2 >= m_phimax ) phi2 -= 2.*M_PI; 

      const int locateid = Locate(r2 , phi2 , z2); //find where the cluster is 
      
      if( z2 > m_zmin || z2 < m_zmax ) GEM_Mod_Arr[locateid-1]++; // the array ath the cluster location - counts the number of clusters in each array, (IFF its in the volume !!!)
      */

      // cluster errors
      const auto cluster_rphi_error = 0.015;
      const auto cluster_z_error = 0.075;
      
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
	  const float r = get_r( projection.x(), projection.y() );
	  const float dr = cluster_r - r;
	  if(h_dca_layer) h_dca_layer->Fill(r, dca);
	  if(h_deltarphi_layer_south && clus_centroid.z() < 0) h_deltarphi_layer_south->Fill(r, drp);
	  if(h_deltarphi_layer_north && clus_centroid.z() > 0) h_deltarphi_layer_north->Fill(r, drp);
	  if(h_deltaz_layer) h_deltaz_layer->Fill(r, dz);
	  if(h_deltar_r) h_deltar_r->Fill(r, dr);
	  if(h_entries)
	    {
	      auto phi = cluster_phi;
	      while( phi < m_phimin ) phi += 2.*M_PI;
	      while( phi >= m_phimax ) phi -= 2.*M_PI;
	      h_entries->Fill( phi, cluster_r, cluster_z );
	    }
	  if(h_xy) h_xy->Fill(clus_centroid.x(), clus_centroid.y());
	  if(h_xz) h_xz->Fill(clus_centroid.x(), clus_centroid.z());
	  if(h_xy_pca) h_xy_pca->Fill(projection.x(), projection.y());
	  if(h_xz_pca) h_xz_pca->Fill(projection.x(), projection.z());
	  if(h_dca_path) h_dca_path->Fill(pathlength,dca);
	  if(h_zr) h_zr->Fill(clus_centroid.z(), cluster_r);
	  if(h_zr_pca) h_zr_pca->Fill(projection.z(), r);
	  if(h_dz_z)h_dz_z->Fill(projection.z(), clus_centroid.z()- projection.z());
          //if(h_clusters)h_clusters->Fill(clus_centroid.x(),clus_centroid.y(),clus_centroid.z());
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
      const auto i = get_cell_index( clus_centroid );
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

// We will eventually bring this back when we fix clustering
/*
  for(int GEMS_iter = 0; GEMS_iter < 72; GEMS_iter++)
  {
    if(GEM_Mod_Arr[GEMS_iter] > 0){  //laser 0
      h_GEMs_hit->Fill(trkid + 0.5);
    }
  }
*/

}


//_____________________________________________________________________
int TpcDirectLaserReconstruction::get_cell_index( const TVector3& global ) const
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
//_____________________________________________________________________
int TpcDirectLaserReconstruction::Locate(float r , float phi , float z)
{

  ////////////////////////////////////////////////////////////////////////
  //          North Label Conventions                                  //
  ///////////////////////////////////////////////////////////////////////
  //   1 - 00_R1   16 - 05_R1   31 - 10_R1    
  //   2 - 00_R2   17 - 05_R2   32 - 10_R2   
  //   3 - 00_R3   18 - 05_R3   33 - 10_R3 
  //   4 - 01_R1   19 - 06_R1   34 - 11_R1
  //   5 - 01_R2   20 - 06_R2   35 - 11_R2
  //   6 - 01_R3   21 - 06_R3   36 - 11_R3
  //   7 - 02_R1   22 - 07_R1
  //   8 - 02_R2   23 - 07_R2
  //   9 - 02_R3   24 - 07_R3
  //  10 - 03_R1   25 - 08_R1
  //  11 - 03_R2   26 - 08_R2
  //  12 - 03_R3   27 - 08_R3
  //  13 - 04_R1   28 - 09_R1
  //  14 - 04_R2   29 - 09_R2
  //  15 - 04_R3   30 - 09_R3

  //////////////////////////////////////////////////////////////////////// 
  //          South Label Conventions                                   //  
  ///////////////////////////////////////////////////////////////////////  
  //  37 - 12_R1   52 - 17_R1   67 - 22_R1    
  //  38 - 12_R2   53 - 17_R2   68 - 22_R2   
  //  39 - 12_R3   54 - 17_R3   69 - 22_R3 
  //  40 - 13_R1   55 - 18_R1   70 - 23_R1
  //  41 - 13_R2   56 - 18_R2   71 - 23_R2
  //  42 - 13_R3   57 - 18_R3   72 - 23_R3
  //  43 - 14_R1   58 - 19_R1
  //  44 - 14_R2   59 - 19_R2
  //  45 - 14_R3   60 - 19_R3
  //  46 - 15_R1   61 - 20_R1
  //  47 - 15_R2   62 - 20_R2
  //  48 - 15_R3   63 - 20_R3
  //  49 - 16_R1   64 - 21_R1
  //  50 - 16_R2   65 - 21_R2
  //  51 - 16_R3   66 - 21_R3


  int GEM_Mod_ID; //integer from 1 to 72

  const float Angle_Bins[13] = {23*M_PI/12,1*M_PI/12,3*M_PI/12,5*M_PI/12,7*M_PI/12,9*M_PI/12,11*M_PI/12,13*M_PI/12,15*M_PI/12,17*M_PI/12,19*M_PI/12,21*M_PI/12,23*M_PI/12}; //12 angle bins on each side

  const float r_bins[4] = {30,46,62,78}; //4 r bins 

  int angle_id = 0; // default to 0
  int r_id = 0; //default to 0
  int side_id = 0; //default to 0 (north side)

  for(int r_iter = 0; r_iter < 3; r_iter++ ){
    if( r > r_bins[r_iter] && r < r_bins[r_iter+1] )
    {
      r_id = r_iter+1;
      break; //break out of the for loop if you found where the hit is in r
    }
  }

  for(int ang_iter = 0; ang_iter < 12; ang_iter++ ){
    if( phi > Angle_Bins[ang_iter] && phi < Angle_Bins[ang_iter+1] ) 
    {
      angle_id = ang_iter;
      break; //break out the for loop if you found where the hit is in phi
    }
  }
  if(z < 0)side_id=1; //(south side)

  return GEM_Mod_ID = (36*side_id) + (3*angle_id) + r_id;
}

float TpcDirectLaserReconstruction::GetRelPhi( float xorig, float yorig, float x, float y, float phiorig ) 
{
// getting the relative phi in the frame rotated to the ith laser origin
//
// inputs: 
//	xorig = x coordinate of ith laser origin
//      yorig = y coordinate of ith laser origin
//          x = x coordinate of arbitrary point in TPC volume
//          y = y coordinate of arbitrary point in TPC volume
//    phiorig = phi (from + x axis) of ith laser origin 
//   
// output:
//     relphi = azimuthal displacement of arbitrary point in TPC volume from ith laser origin

  if (phiorig < 0) phiorig += 2. * M_PI;
  else if (phiorig > 2.*M_PI) phiorig -= 2. * M_PI;

  float dx = x - xorig;
  float dy = y - yorig;
  float relphi = atan2(dy, dx) - phiorig;
  if (relphi < 0) relphi += 2. * M_PI;
  else if (relphi > 2.*M_PI) relphi -= 2. * M_PI;

  return relphi;
}

float TpcDirectLaserReconstruction::GetRelTheta( float xorig, float yorig, float zorig, float x, float y, float z, float thetaorig, float phiorig )
{
// getting the relative theta in the frame rotated to the ith laser origin
//
// inputs: 
//	xorig = x coordinate of ith laser origin
//      yorig = y coordinate of ith laser origin
//      zorig = z cooridante of ith laser origin
//          x = x coordinate of arbitrary point in TPC volume
//          y = y coordinate of arbitrary point in TPC volume
//          z = z coordinate of arbitrary point in TPC volume
//  thetaorig = theta (from + y axis?) of ithe laser origin
//    phiorig = phi (from + x axis) of ith laser origin 
//   
// output:
//     reltheta = polar displacement of arbitrary point in TPC volume from ith laser origin

  if (phiorig < 0) phiorig += 2. * M_PI;
  else if (phiorig > 2.*M_PI) phiorig -= 2. * M_PI;

  if (thetaorig < 0) thetaorig += 2. * M_PI;
  else if (thetaorig > 2.*M_PI) thetaorig -= 2. * M_PI;


  float dx = x - xorig;
  float dy = y - yorig;
  float dz = z - zorig;
  float r = sqrt(dx*dx + dy*dy + dz*dz);

  float cos_beta = (dx*cos(phiorig) + dy*sin(phiorig))/r;
  float sin_beta = dz/r;

  float reltheta = acos(cos_beta*cos(thetaorig) + sin_beta*sin(thetaorig)) - M_PI/2.;
  if (reltheta < 0) reltheta += 2. * M_PI;
  else if (reltheta > 2.*M_PI) reltheta -= 2. * M_PI;

  return reltheta;
}

bool TpcDirectLaserReconstruction::sameSign(float num1, float num2)
{
    return (num1 >= 0 && num2 >= 0) || (num1 < 0 && num2 < 0);
}


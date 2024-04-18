/**
 * \file TpcSpaceChargeMatrixInversion.cc
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeMatrixInversion.h"
#include "TpcSpaceChargeReconstructionHelper.h"
#include "TpcSpaceChargeMatrixContainerv1.h"

#include <frog/FROG.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <TFile.h>
#include <TH3.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <memory>

namespace
{
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
TpcSpaceChargeMatrixInversion::TpcSpaceChargeMatrixInversion( const std::string& name ):
  Fun4AllBase( name)
{}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion::add_from_file( const std::string& shortfilename, const std::string& objectname )
{
  // get filename from frog
  FROG frog;
  const auto filename = frog.location( shortfilename );

  // open TFile
  std::unique_ptr<TFile> inputfile( TFile::Open( filename ) );
  if( !inputfile )
  {
    std::cout << "TpcSpaceChargeMatrixInversion::add_from_file - could not open file " << filename << std::endl;
    return false;
  }

  // load object from input file
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source( dynamic_cast<TpcSpaceChargeMatrixContainer*>( inputfile->Get( objectname.c_str() ) ) );
  if( !source )
  {
    std::cout << "TpcSpaceChargeMatrixInversion::add_from_file - could not find object name " << objectname << " in file " << filename << std::endl;
    return false;
  }

  // add object
  return add( *source.get() );
}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion::add( const TpcSpaceChargeMatrixContainer& source )
{
  // check internal container, create if necessary
  if( !m_matrix_container )
  {
    m_matrix_container.reset( new TpcSpaceChargeMatrixContainerv1 );

    // get grid dimensions from source
    int phibins = 0;
    int rbins = 0;
    int zbins = 0;
    source.get_grid_dimensions( phibins, rbins, zbins );

    // assign
    m_matrix_container->set_grid_dimensions( phibins, rbins, zbins );
  }

  // add content
  return m_matrix_container->add( source );
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion::calculate_distortions()
{

  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );

  // create output histograms
  std::unique_ptr<TH3> hentries( new TH3F( "hentries_rec", "hentries_rec", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax ) );
  std::unique_ptr<TH3> hphi( new TH3F( "hDistortionP_rec", "hDistortionP_rec", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax ) );
  std::unique_ptr<TH3> hz( new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax ) );
  std::unique_ptr<TH3> hr( new TH3F( "hDistortionR_rec", "hDistortionR_rec", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax ) );

  // set axis labels
  for( const auto& h:{ hentries.get(), hphi.get(), hz.get(), hr.get() } )
  {
    h->GetXaxis()->SetTitle( "#phi (rad)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "z (cm)" );
  }

  // matrix convenience definition
  /* number of coordinates must match that of the matrix container */
  static constexpr int ncoord = 3;
  using matrix_t = Eigen::Matrix<float, ncoord, ncoord >;
  using column_t = Eigen::Matrix<float, ncoord, 1 >;

  // loop over bins
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {

    // get cell index
    const auto icell = m_matrix_container->get_cell_index( iphi, ir, iz );

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 2;
    const auto cell_entries = m_matrix_container->get_entries(icell);
    if( cell_entries < min_cluster_count ) continue;

    // build eigen matrices from container
    matrix_t lhs;
    for( int i = 0; i < ncoord; ++i )
      for( int j = 0; j < ncoord; ++j )
    { lhs(i,j) = m_matrix_container->get_lhs( icell, i, j ); }

    column_t rhs;
    for( int i = 0; i < ncoord; ++i )
    { rhs(i) = m_matrix_container->get_rhs( icell, i ); }

    if (Verbosity())
    {
      // print matrices and entries
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - inverting bin " << iz << ", " << ir << ", " << iphi << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - entries: " << cell_entries << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - lhs: \n" << lhs << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - rhs: \n" << rhs << std::endl;
    }

    // calculate result using linear solving
    const auto cov = lhs.inverse();
    auto partialLu = lhs.partialPivLu();
    const auto result = partialLu.solve( rhs );

    // fill histograms
    hentries->SetBinContent( iphi+1, ir+1, iz+1, cell_entries );

    hphi->SetBinContent( iphi+1, ir+1, iz+1, result(0) );
    hphi->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(0,0) ) );

    hz->SetBinContent( iphi+1, ir+1, iz+1, result(1) );
    hz->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(1,1) ) );

    hr->SetBinContent( iphi+1, ir+1, iz+1, result(2) );
    hr->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(2,2) ) );

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - drphi: " << result(0) << " +/- " << std::sqrt( cov(0,0) ) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - dz: " << result(1) << " +/- " << std::sqrt( cov(1,1) ) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion::calculate_distortions - dr: " << result(2) << " +/- " << std::sqrt( cov(2,2) ) << std::endl;
      std::cout << std::endl;
    }
  }

  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = []( TH3* h, const TString& name )
  {
    const auto& result = TpcSpaceChargeReconstructionHelper::split( h );
    std::unique_ptr<TH3> hneg( std::get<0>(result) );
    std::unique_ptr<TH3> hpos( std::get<1>(result) );

    return std::make_tuple(
      TpcSpaceChargeReconstructionHelper::add_guarding_bins( hneg.get(), Form( "%s_negz", name.Data() ) ),
      TpcSpaceChargeReconstructionHelper::add_guarding_bins( hpos.get(), Form( "%s_posz", name.Data() ) ));
  };


  // apply finishing transformations to histograms and save in container
  if( !m_dcc_average ) { m_dcc_average.reset(new TpcDistortionCorrectionContainer); }

  std::tie( m_dcc_average->m_hentries[0], m_dcc_average->m_hentries[1] ) = process_histogram( hentries.get(), "hentries" );
  std::tie( m_dcc_average->m_hDRint[0], m_dcc_average->m_hDRint[1] ) = process_histogram( hr.get(), "hIntDistortionR" );
  std::tie( m_dcc_average->m_hDPint[0], m_dcc_average->m_hDPint[1] ) = process_histogram( hphi.get(), "hIntDistortionP" );
  std::tie( m_dcc_average->m_hDZint[0], m_dcc_average->m_hDZint[1] ) = process_histogram( hz.get(), "hIntDistortionZ" );

}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion::load_cm_distortion_correction( const std::string& filename )
{
  // open TFile
  auto distortion_tfile = TFile::Open(filename.c_str());
  if( !distortion_tfile )
  {
    std::cout << "TpcSpaceChargeMatrixInversion::load_cm_distortion_correction - cannot open " << filename << std::endl;
    exit(1);
  }

  // make sure container exists
  if( !m_dcc_cm ) { m_dcc_cm.reset(new TpcDistortionCorrectionContainer); }

  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
  for (int j = 0; j < 2; ++j)
  {
    m_dcc_cm->m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
    assert(m_dcc_cm->m_hDPint[j]);
    m_dcc_cm->m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
    assert(m_dcc_cm->m_hDRint[j]);
    m_dcc_cm->m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
    assert(m_dcc_cm->m_hDZint[j]);
  }

  // check histogram dimension. We expect 2D histograms
  m_dcc_cm->m_dimensions = m_dcc_cm->m_hDPint[0]->GetDimension();
  assert( m_dcc_cm->m_dimensions == 2 );
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion::save_distortions(const std::string& filename)
{
  if( !m_dcc_average )
  {
    std::cout << "TpcSpaceChargeMatrixInversion::save_distortions - invalid histograms." << std::endl;
    return;
  }

  // save everything to root file
  std::cout << "TpcSpaceChargeMatrixInversion::save_distortions - writing histograms to " << filename << std::endl;
  std::unique_ptr<TFile> outputfile( TFile::Open( filename.c_str(), "RECREATE" ) );
  outputfile->cd();

  for( const auto& h: {m_dcc_average->m_hentries, m_dcc_average->m_hDRint, m_dcc_average->m_hDPint, m_dcc_average->m_hDZint } )
  {
    if( h[0] ) h[0]->Write();
    if( h[1] ) h[1]->Write();
  }

  // close TFile
  outputfile->Close();
}

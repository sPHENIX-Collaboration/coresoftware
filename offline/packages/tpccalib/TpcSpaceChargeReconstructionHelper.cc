/**
 * \file TpcSpaceChargeReconstructionHelper.cc
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeReconstructionHelper.h"

#include <micromegas/CylinderGeomMicromegas.h>

#include <TH3.h>
#include <TString.h>

#include <iostream>

namespace
{
  /// square
  template<class T> 
  inline constexpr T square( T x ) { return x*x; }

  // regularize angle between 0 and 2PI
  template<class T>
  inline constexpr T get_bound_angle( T phi ) 
  {
    while( phi < 0 ) phi += 2*M_PI;
    while( phi >= 2*M_PI ) phi -= 2*M_PI;
    return phi;
  }
 
  // angle from a given sector
  constexpr double get_sector_phi( int isec ) { return isec*M_PI/6; }

  // Micromegas geometry
  // TODO: should get those numbers from actual geometry configuration

  /// fully equiped sector
  static constexpr double isec_ref = 9;
  static constexpr double phi_ref = get_sector_phi(isec_ref);

  // radius of the outermost micromegas layer
  static constexpr double r_ref = 85.1;

  /// micromegas azimutal angle. It is used for interpolation between sectors
  /** 
   * 31.6cm corresponds to the micromegas tile width from CAD drawings, also used in PHG4MicromegasDetector.cc
   * there is a reduction factor of 0.6, to avoid side effects
   */
  static constexpr double delta_phi_mm = 0.6*(31.6/CylinderGeomMicromegas::reference_radius); 
  
  /// z extrapolation window
  static constexpr double zextrap_min = 53.2 - 5.0;
  static constexpr double zextrap_max = 56.6 + 5.0;
  
}

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_z( TH3* hin )
{
  if( !hin ) return;

  // get reference phi bin
  const int phibin_ref = hin->GetXaxis()->FindBin( phi_ref );

  // loop over radial bins
  for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
  {

    // get current radius
    const auto r = hin->GetYaxis()->GetBinCenter( ir+1 );

    // get z integration window for reference
    const auto zextrap_min_loc = zextrap_min * r/r_ref;
    const auto zextrap_max_loc = zextrap_max * r/r_ref;

    // get corresponding bins
    const int zbin_min[2] = { hin->GetZaxis()->FindBin( -zextrap_max_loc ), hin->GetZaxis()->FindBin( zextrap_min_loc ) };
    const int zbin_max[2] = { hin->GetZaxis()->FindBin( -zextrap_min_loc ), hin->GetZaxis()->FindBin( zextrap_max_loc ) };

    for( int isign = 0; isign < 2; ++isign )
    {
      // adjust z positions
      const auto z_min = hin->GetZaxis()->GetBinCenter( zbin_min[isign] );
      const auto z_max = hin->GetZaxis()->GetBinCenter( zbin_max[isign] );

      // get reference
      const auto content_min = hin->GetBinContent( phibin_ref, ir+1, zbin_min[isign] );
      const auto content_max = hin->GetBinContent( phibin_ref, ir+1, zbin_max[isign] );
      const auto error_min = hin->GetBinError( phibin_ref, ir+1, zbin_min[isign] );
      const auto error_max = hin->GetBinError( phibin_ref, ir+1, zbin_max[isign] );

      // loop over z bins
      for( int iz = zbin_min[isign]+1; iz < zbin_max[isign]; ++iz )
      {

        const auto z = hin->GetZaxis()->GetBinCenter( iz );

        // interpolate
        const auto alpha_min = (z_max-z)/(z_max-z_min);
        const auto alpha_max = (z-z_min)/(z_max-z_min);

        const auto content = alpha_min*content_min + alpha_max*content_max;
        const auto error = std::sqrt(square( alpha_min * error_min ) + square( alpha_max*error_max));

        hin->SetBinContent( phibin_ref, ir+1, iz, content );
        hin->SetBinError( phibin_ref, ir+1, iz, error );
      }
    }
  }
}

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_phi1( TH3* hin )
{
  if( !hin ) return;
  
  // get phi bin range for a given sector
  auto get_phibin_range = []( TH3* hin, int isec  )
  {

    // get corresponding first and last bin
    const double phi = get_sector_phi( isec );
    const double phi_min = get_bound_angle( phi - M_PI/12 );
    const double phi_max = get_bound_angle( phi + M_PI/12 ); 

    // find corresponding bins
    const int phibin_min = hin->GetXaxis()->FindBin( phi_min );
    const int phibin_max = hin->GetXaxis()->FindBin( phi_max );    
    return std::make_pair( phibin_min, phibin_max ); 
  };
  
  // get reference bins
  const auto [phibin_min_ref, phibin_max_ref] = get_phibin_range( hin, isec_ref );
  
  // get number of phi bins
  const int nphibins = hin->GetNbinsX();

  // copy all r and z bins from reference phi bin to destination
  auto copy_phi_bin = []( TH3* hin, int phibin_ref, int phibin_dest )
  {
    
    // loop over radial bins
    for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
    {
      
      // loop over z bins
      for( int iz = 0; iz < hin->GetZaxis()->GetNbins(); ++iz )
      {
        const auto content_ref = hin->GetBinContent( phibin_ref, ir+1, iz+1 );
        const auto error_ref = hin->GetBinError( phibin_ref, ir+1, iz+1 );
        
        // calculate scale factor
        const auto scale = 1;
        
        // assign to output histogram
        hin->SetBinContent( phibin_dest, ir+1, iz+1, content_ref*scale );
        hin->SetBinError( phibin_dest, ir+1, iz+1, error_ref*std::abs(scale) );
      }
    }
    
  };
  
  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {
    // skip reference sector
    if( isec == isec_ref ) continue;

    const auto [phibin_min, phibin_max] = get_phibin_range( hin, isec );

    // loop over bins
    for( int ibin = 0; ibin < phibin_max_ref - phibin_min_ref; ++ibin )
    {
      int phibin_ref = (phibin_min_ref + ibin); 
      if( phibin_ref > nphibins ) phibin_ref -= nphibins;
      
      int phibin = (phibin_min + ibin);   
      if( phibin > nphibins ) phibin -= nphibins;
      
      copy_phi_bin( hin, phibin_ref, phibin );
    }
  }
}

//_______________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_phi2( TH3* hin )
{
  if( !hin ) return;

  
  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {    
    // get phi range for interpolation from this sector to the next
    const double phi_min = get_sector_phi(isec)+delta_phi_mm/2;
    const double phi_max = get_sector_phi(isec+1)-delta_phi_mm/2;
    
    // get corresponding bins
    const int phibin_min = hin->GetXaxis()->FindBin(get_bound_angle(phi_min));
    const int phibin_max = hin->GetXaxis()->FindBin(get_bound_angle(phi_max));
        
    // loop over radial bins
    for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
    {
      
      // loop over z bins
      for( int iz = 0; iz < hin->GetZaxis()->GetNbins(); ++iz )
      {
        const auto content_min = hin->GetBinContent( phibin_min, ir+1, iz+1 );
        const auto content_max = hin->GetBinContent( phibin_max, ir+1, iz+1 );
        const auto error_min = hin->GetBinError( phibin_min, ir+1, iz+1 );
        const auto error_max = hin->GetBinError( phibin_max, ir+1, iz+1 );

        // loop over relevant phi bins
        for( int iphi = phibin_min+1; iphi < phibin_max; ++iphi )
        {
          // get phi
          const auto phi = hin->GetXaxis()->GetBinCenter( iphi );
          
          // perform linear extrapolation
          const auto alpha_min = (phi_max-phi)/(phi_max-phi_min);
          const auto alpha_max = (phi-phi_min)/(phi_max-phi_min);
          
          const auto content = alpha_min*content_min + alpha_max*content_max;
          const auto error = std::sqrt(square( alpha_min * error_min ) + square( alpha_max*error_max));
          
          hin->SetBinContent( iphi, ir+1, iz+1, content );
          hin->SetBinError( iphi, ir+1, iz+1, error );
        }
      }
    }
  }
  
}

//_______________________________________________
std::tuple<TH3*, TH3*> TpcSpaceChargeReconstructionHelper::split( TH3* hin )
{
  if( !hin ) return std::make_tuple<TH3*, TH3*>( nullptr, nullptr );

  auto xaxis = hin->GetXaxis();
  auto yaxis = hin->GetYaxis();
  auto zaxis = hin->GetZaxis();
  auto ibin = zaxis->FindBin( (double) 0 );

  // create histograms
  auto hneg = new TH3F(
    Form( "%s_negz", hin->GetName() ), Form( "%s_negz", hin->GetTitle() ),
    xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax(),
    yaxis->GetNbins(), yaxis->GetXmin(), yaxis->GetXmax(),
    ibin-1, zaxis->GetXmin(), zaxis->GetBinUpEdge( ibin-1 ) );

  auto hpos = new TH3F(
    Form( "%s_posz", hin->GetName() ), Form( "%s_posz", hin->GetTitle() ),
    xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax(),
    yaxis->GetNbins(), yaxis->GetXmin(), yaxis->GetXmax(),
    zaxis->GetNbins() - (ibin-1), zaxis->GetBinLowEdge(ibin), zaxis->GetXmax() );

  // copy content and errors
  for( int ix = 0; ix < xaxis->GetNbins(); ++ix )
    for( int iy = 0; iy < yaxis->GetNbins(); ++iy )
    for( int iz = 0; iz < zaxis->GetNbins(); ++iz )
  {
    const auto content = hin->GetBinContent( ix+1, iy+1, iz+1 );
    const auto error = hin->GetBinError( ix+1, iy+1, iz+1 );

    if( iz < ibin-1 )
    {
      hneg->SetBinContent( ix+1, iy+1, iz+1, content );
      hneg->SetBinError( ix+1, iy+1, iz+1, error );
    } else {
      hpos->SetBinContent( ix+1, iy+1, iz - (ibin-1) + 1, content );
      hpos->SetBinError( ix+1, iy+1, iz - (ibin-1) + 1, error );
    }
  }

  // also copy axis titles
  for( const auto h: {hneg, hpos} )
  {
    h->GetXaxis()->SetTitle( hin->GetXaxis()->GetTitle() );
    h->GetYaxis()->SetTitle( hin->GetYaxis()->GetTitle() );
    h->GetZaxis()->SetTitle( hin->GetZaxis()->GetTitle() );
  }

  return std::make_tuple( hneg, hpos );
}

//___________________________________________________________________________
TH3* TpcSpaceChargeReconstructionHelper::copy_histogram( TH3* hin, const TString& name )
{
  std::array<int, 3> bins;
  std::array<double, 3> x_min;
  std::array<double, 3> x_max;

  int index = 0;
  for( const auto axis:{ hin->GetXaxis(), hin->GetYaxis(), hin->GetZaxis() } )
  {
    // calculate bin width
    const auto bin_width = (axis->GetXmax() - axis->GetXmin())/axis->GetNbins();

    // increase the number of bins by two
    bins[index] = axis->GetNbins()+2;

    // update axis limits accordingly
    x_min[index] = axis->GetXmin()-bin_width;
    x_max[index] = axis->GetXmax()+bin_width;
    ++index;
  }

  // create new histogram
  auto hout = new TH3F( name, name,
    bins[0], x_min[0], x_max[0],
    bins[1], x_min[1], x_max[1],
    bins[2], x_min[2], x_max[2] );

  // update axis legend
  hout->GetXaxis()->SetTitle( hin->GetXaxis()->GetTitle() );
  hout->GetYaxis()->SetTitle( hin->GetYaxis()->GetTitle() );
  hout->GetZaxis()->SetTitle( hin->GetZaxis()->GetTitle() );

  // copy content
  const auto phibins = hin->GetXaxis()->GetNbins();
  const auto rbins = hin->GetYaxis()->GetNbins();
  const auto zbins = hin->GetZaxis()->GetNbins();

  // fill center
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {
    hout->SetBinContent( iphi+2, ir+2, iz+2, hin->GetBinContent( iphi+1, ir+1, iz+1 ) );
    hout->SetBinError( iphi+2, ir+2, iz+2, hin->GetBinError( iphi+1, ir+1, iz+1 ) );
  }

  // fill guarding phi bins
  /*
   * we use 2pi periodicity to do that:
   * - last valid bin is copied to first guarding bin;
   * - first valid bin is copied to last guarding bin
   */
  for( int ir = 0; ir < rbins+2; ++ir )
    for( int iz = 0; iz < zbins+2; ++iz )
  {
    // copy last bin to first guarding bin
    hout->SetBinContent( 1, ir+1, iz+1, hout->GetBinContent( phibins+1, ir+1, iz+1 ) );
    hout->SetBinError( 1, ir+1, iz+1, hout->GetBinError( phibins+1, ir+1, iz+1 ) );

    // copy first bin to last guarding bin
    hout->SetBinContent( phibins+2, ir+1, iz+1, hout->GetBinContent( 2, ir+1, iz+1 ) );
    hout->SetBinError( phibins+2, ir+1, iz+1, hout->GetBinError( 2, ir+1, iz+1 ) );
  }

  // fill guarding r bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int iz = 0; iz < zbins+2; ++iz )
  {
    hout->SetBinContent( iphi+1, 1, iz+1, hout->GetBinContent( iphi+1, 2, iz+1 ) );
    hout->SetBinError( iphi+1, 1, iz+1, hout->GetBinError( iphi+1, 2, iz+1 ) );

    hout->SetBinContent( iphi+1, rbins+2, iz+1, hout->GetBinContent( iphi+1, rbins+1, iz+1 ) );
    hout->SetBinError( iphi+1, rbins+2, iz+1, hout->GetBinError( iphi+1, rbins+1, iz+1 ) );
  }

  // fill guarding z bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int ir = 0; ir < rbins+2; ++ir )
  {
    hout->SetBinContent( iphi+1, ir+1, 1, hout->GetBinContent( iphi+1, ir+1, 2 ) );
    hout->SetBinError( iphi+1, ir+1, 1, hout->GetBinError( iphi+1, ir+1, 2 ) );

    hout->SetBinContent( iphi+1, ir+1, zbins+2, hout->GetBinContent( iphi+1, ir+1, zbins+1 ) );
    hout->SetBinError( iphi+1, ir+1, zbins+2, hout->GetBinError( iphi+1, ir+1, zbins+1 ) );
  }

  return hout;

}

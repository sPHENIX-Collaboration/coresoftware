#include "TpcSpaceChargeReconstructionHelper.h"

#include <TH3.h>
#include <TString.h>

namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// Micromegas geometry
  /// TODO: should get those numbers from actual geometry configuration
  // fully equiped sector
  static constexpr double isec_ref = 3;
  static constexpr double phi_ref = isec_ref*M_PI/6 + M_PI/12;

  // radius of the micromegas layer
  static constexpr double r_ref = 82;

  // z extrapolation window
  static constexpr double zextrap_min = 48;
  static constexpr double zextrap_max = 58;

  // Micromegas acceptance in incomplete sectors
  static constexpr double zref = 33.25;
  static constexpr double length = 50 - 5;
  static constexpr double zref_min = zref - length/2;
  static constexpr double zref_max = zref + length/2;

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

  // get reference phi bin
  const int phibin_ref = hin->GetXaxis()->FindBin( phi_ref );

  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {

    // skip reference sector
    if( isec == isec_ref ) continue;

    // get relevant phi and corresponding bin
    const double phi = isec*M_PI/6 + M_PI/12;
    const int phibin = hin->GetXaxis()->FindBin( phi );

    // loop over radial bins
    for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
    {

      // get current radius
      const auto r = hin->GetYaxis()->GetBinCenter( ir+1 );

      // get z integration window for reference
      const auto zref_min_loc = zref_min * r/r_ref;
      const auto zref_max_loc = zref_max * r/r_ref;

      // get corresponding bins
      const int zbin_ref_neg[2] = { hin->GetZaxis()->FindBin( -zref_max_loc ), hin->GetZaxis()->FindBin( -zref_min_loc ) };
      const int zbin_ref_pos[2] = { hin->GetZaxis()->FindBin( zref_min_loc ), hin->GetZaxis()->FindBin( zref_max_loc ) };

      // loop over z bins
      for( int iz = 0; iz < hin->GetZaxis()->GetNbins(); ++iz )
      {
        const auto content_ref = hin->GetBinContent( phibin_ref, ir+1, iz+1 );
        const auto error_ref = hin->GetBinError( phibin_ref, ir+1, iz+1 );

        #if true
        // calculate scale factor
        const auto z = hin->GetZaxis()->GetBinCenter( iz+1 );
        const auto norm_ref = hin->Integral( phibin_ref, phibin_ref, ir+1, ir+1, (z>0) ? zbin_ref_pos[0]:zbin_ref_neg[0], (z>0) ? zbin_ref_pos[1]:zbin_ref_neg[1] );
        const auto norm_loc = hin->Integral( phibin, phibin, ir+1, ir+1, (z>0) ? zbin_ref_pos[0]:zbin_ref_neg[0], (z>0) ? zbin_ref_pos[1]:zbin_ref_neg[1] );
        const auto scale = (norm_ref == 0) ? 1:norm_loc/norm_ref;

        #else
        const auto scale = 1;
        #endif

        // assign to output histogram
        hin->SetBinContent( phibin, ir+1, iz+1, content_ref*scale );
        hin->SetBinError( phibin, ir+1, iz+1, error_ref*scale );
      }
    }
  }
}

//_______________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_phi2( TH3* hin )
{
  if( !hin ) return;

  for( int iphi = 0; iphi < hin->GetXaxis()->GetNbins(); ++iphi )
  {

    // find nearest sector phi bins
    const auto phi = hin->GetXaxis()->GetBinCenter( iphi+1 );
    const int isec = std::floor( (phi - M_PI/12)/(M_PI/6) );
    double phi_min =  isec*M_PI/6 + M_PI/12;
    double phi_max =  phi_min + M_PI/6;

    if( phi_min < 0 ) phi_min += 2*M_PI;
    if( phi_max >= 2*M_PI ) phi_max -= 2*M_PI;

    const auto phibin_min = hin->GetXaxis()->FindBin( phi_min );
    if( phibin_min == iphi+1 ) continue;

    const auto phibin_max = hin->GetXaxis()->FindBin( phi_max );
    if( phibin_max == iphi+1 ) continue;

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

        // perform linear extrapolation
        const auto alpha_min = (phi_max-phi)/(phi_max-phi_min);
        const auto alpha_max = (phi-phi_min)/(phi_max-phi_min);

        const auto content = alpha_min*content_min + alpha_max*content_max;
        const auto error = std::sqrt(square( alpha_min * error_min ) + square( alpha_max*error_max));

        hin->SetBinContent( iphi+1, ir+1, iz+1, content );
        hin->SetBinError( iphi+1, ir+1, iz+1, error );
      }

    }
  }
}

//_______________________________________________
std::array<TH3*, 2> TpcSpaceChargeReconstructionHelper::split( TH3* hin )
{
  if( !hin ) return {{nullptr, nullptr}};

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

  return {{hneg, hpos}};
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
    hout->SetBinError( iphi+1, rbins+1, iz+1, hout->GetBinError( iphi+1, rbins+1, iz+1 ) );
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

/**
 * \file TpcSpaceChargeReconstructionHelper.cc
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeReconstructionHelper.h"

#include <TH2.h>
#include <TH3.h>
#include <TString.h>

#include <iostream>

namespace
{
  /// square
  template <class T>
  inline constexpr T square(T x)
  {
    return x * x;
  }

  // regularize angle between 0 and 2PI
  template <class T>
  inline constexpr T get_bound_angle(T phi)
  {
    while (phi < 0)
    {
      phi += 2 * M_PI;
    }
    while (phi >= 2 * M_PI)
    {
      phi -= 2 * M_PI;
    }
    return phi;
  }

  /// shortcut to angular window, needed to define TPOT acceptance
  using range_t = std::pair<double, double>;

  /// list of angular windows
  using range_list_t = std::vector<range_t>;

  /// make sure angles in a given window are betwwen [0, 2PI]
  range_t transform_range(const range_t& a)
  {
    return {get_bound_angle(a.first), get_bound_angle(a.second)};
  }

  ///@name detector geometry
  //@{
  /// phi range for central sector (TPC sectors 9 and 21)
  const range_t phi_range_central = transform_range({-1.742, -1.43979});

  /// phi range for east sector (TPC sectors 8 and 20)
  const range_t phi_range_east = transform_range({-2.27002, -1.9673});

  /// phi range for west sector (TPC sectors 10 and 22)
  const range_t phi_range_west = transform_range({-1.21452, -0.911172});

  /// list of theta (polar) angles for each micromeas in central sector
  const range_list_t theta_range_central = {{-0.918257, -0.613136}, {-0.567549, -0.031022}, {0.0332154, 0.570419}, {0.613631, 0.919122}};

  /// list of theta (polar) angles for each micromeas in east sector
  const range_list_t theta_range_east = {{-0.636926, -0.133603}, {0.140678, 0.642714}};

  /// list of theta (polar) angles for each micromeas in west sector
  const range_list_t theta_range_west = {{-0.643676, -0.141004}, {0.13485, 0.640695}};
  //@}

  /// short class to check if a given value is in provided range
  class range_ftor_t
  {
   public:
    explicit range_ftor_t(double value)
      : m_value(value){};
    bool operator()(const range_t& range)
    {
      return m_value > range.first && m_value < range.second;
    }

   private:
    double m_value = 0;
  };

  /// returns true if given value is in provided range
  bool in_range(const double& value, const range_t& range)
  {
    return range_ftor_t(value)(range);
  }

  /// returns true if given value is in any of the provided range
  bool in_range(const double& value, const range_list_t& range_list)
  {
    return std::any_of(range_list.begin(), range_list.end(), range_ftor_t(value));
  }

}  // namespace

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::create_tpot_mask(TH3* hmask)
{
  hmask->Reset();

  // loop over bins
  for (int ip = 0; ip < hmask->GetNbinsX(); ++ip)
  {
    for (int ir = 0; ir < hmask->GetNbinsY(); ++ir)
    {
      for (int iz = 0; iz < hmask->GetNbinsZ(); ++iz)
      {
        const double phi = hmask->GetXaxis()->GetBinCenter(ip + 1);
        const double r = hmask->GetYaxis()->GetBinCenter(ir + 1);
        const double z = hmask->GetZaxis()->GetBinCenter(iz + 1);
        const double theta = std::atan2(z, r);
        if (in_range(phi, phi_range_central) && in_range(theta, theta_range_central))
        {
          hmask->SetBinContent(ip + 1, ir + 1, iz + 1, 1);
        }
      }
    }
  }
}

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_z(TH3* source, const TH3* mask)
{
  if (!source)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_z - invalid source histogram" << std::endl;
    return;
  }

  if (!mask)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_z - invalid mask histogram" << std::endl;
    return;
  }

  // loop over phi bins
  for (int ir = 0; ir < source->GetNbinsY(); ++ir)
  {
    for (int ip = 0; ip < source->GetNbinsX(); ++ip)
    {
      int iz_min = -1;
      for (int iz = 0; iz < source->GetNbinsZ(); ++iz)
      {
        bool in_range = mask->GetBinContent(ip + 1, ir + 1, iz + 1) > 0;
        if (in_range)
        {
          iz_min = iz;
          continue;
        }

        // iz is not in range. Check if iz_min was set. If not, skip this bin. If yes, find next bin in range
        if (iz_min < 0)
        {
          continue;
        }

        int iz_max = iz + 1;
        for (; iz_max < source->GetNbinsZ(); ++iz_max)
        {
          in_range = mask->GetBinContent(ip + 1, ir + 1, iz_max + 1) > 0;
          if (in_range)
          {
            break;
          }
        }

        // check interpolation range validity
        if (iz_max == source->GetNbinsZ())
        {
          continue;
        }

        // do the interpolation
        const double z_min = source->GetZaxis()->GetBinCenter(iz_min + 1);
        const double z_max = source->GetZaxis()->GetBinCenter(iz_max + 1);
        const double z = source->GetZaxis()->GetBinCenter(iz + 1);
        const double alpha_min = (z_max - z) / (z_max - z_min);
        const double alpha_max = (z - z_min) / (z_max - z_min);

        const double distortion_min = source->GetBinContent(ip + 1, ir + 1, iz_min + 1);
        const double distortion_max = source->GetBinContent(ip + 1, ir + 1, iz_max + 1);
        const double distortion = alpha_min * distortion_min + alpha_max * distortion_max;

        const double error_min = source->GetBinError(ip + 1, ir + 1, iz_min + 1);
        const double error_max = source->GetBinError(ip + 1, ir + 1, iz_max + 1);
        const double error = std::sqrt(square(alpha_min * error_min) + square(alpha_max * error_max));

        source->SetBinContent(ip + 1, ir + 1, iz + 1, distortion);
        source->SetBinError(ip + 1, ir + 1, iz + 1, error);
      }
    }
  }
}

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_z2(TH3* source, const TH3* mask, uint8_t side)
{
  if (!source)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_z - invalid source histogram" << std::endl;
    return;
  }

  if (!mask)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_z - invalid mask histogram" << std::endl;
    return;
  }

  // loop over phi bins
  for (int ir = 0; ir < source->GetNbinsY(); ++ir)
  {
    for (int ip = 0; ip < source->GetNbinsX(); ++ip)
    {
      // handle pad plane on the positive side (zmax)
      if (side & Side_positive)
      {
        // find first non empty bin, starting from last zaxis bin
        int iz_min = -1;
        for (int iz = source->GetNbinsZ() - 1; iz >= 0; --iz)
        {
          const bool in_range = mask->GetBinContent(ip + 1, ir + 1, iz + 1) > 0;
          if (in_range)
          {
            iz_min = iz;
            break;
          }
        }

        // check bin validity
        if (iz_min < 0)
        {
          continue;
        }

        // loop over empty bins and do the interpolation
        for (int iz = iz_min + 1; iz < source->GetNbinsZ(); ++iz)
        {
          const double z_min = source->GetZaxis()->GetBinCenter(iz_min + 1);
          const double z_max = source->GetZaxis()->GetXmax();
          const double z = source->GetZaxis()->GetBinCenter(iz + 1);
          const double alpha_min = (z_max - z) / (z_max - z_min);
          const double alpha_max = (z - z_min) / (z_max - z_min);

          const double distortion_min = source->GetBinContent(ip + 1, ir + 1, iz_min + 1);
          const double distortion_max = 0;
          const double distortion = alpha_min * distortion_min + alpha_max * distortion_max;

          const double error_min = source->GetBinError(ip + 1, ir + 1, iz_min + 1);
          const double error_max = 0;
          const double error = std::sqrt(square(alpha_min * error_min) + square(alpha_max * error_max));

          source->SetBinContent(ip + 1, ir + 1, iz + 1, distortion);
          source->SetBinError(ip + 1, ir + 1, iz + 1, error);
        }
      }

      // handle pad plane on the negative side (zmin
      if (side & Side_negative)
      {
        // find first non empty bin, starting from first zaxis bin
        int iz_max = -1;
        for (int iz = 0; iz < source->GetNbinsZ(); ++iz)
        {
          const bool in_range = mask->GetBinContent(ip + 1, ir + 1, iz + 1) > 0;
          if (in_range)
          {
            iz_max = iz;
            break;
          }
        }

        // check bin validity
        if (iz_max < 0)
        {
          continue;
        }

        // loop over empty bins and do the interpolation
        for (int iz = 0; iz < iz_max; ++iz)
        {
          const double z_min = source->GetZaxis()->GetXmin();
          const double z_max = source->GetZaxis()->GetBinCenter(iz_max + 1);
          const double z = source->GetZaxis()->GetBinCenter(iz + 1);
          const double alpha_min = (z_max - z) / (z_max - z_min);
          const double alpha_max = (z - z_min) / (z_max - z_min);

          const double distortion_min = 0;
          const double distortion_max = source->GetBinContent(ip + 1, ir + 1, iz_max + 1);
          const double distortion = alpha_min * distortion_min + alpha_max * distortion_max;

          const double error_min = 0;
          const double error_max = source->GetBinError(ip + 1, ir + 1, iz_max + 1);
          const double error = std::sqrt(square(alpha_min * error_min) + square(alpha_max * error_max));

          source->SetBinContent(ip + 1, ir + 1, iz + 1, distortion);
          source->SetBinError(ip + 1, ir + 1, iz + 1, error);
        }
      }
    }
  }
}

//____________________________________________________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_phi1(TH3* source, const TH2* source_cm, const TH3* mask)
{
  if (!source)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_phi1 - invalid source histogram" << std::endl;
    return;
  }

  if (!mask)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_phi1 - invalid mask histogram" << std::endl;
    return;
  }

  // loop over phi bins
  for (int ip = 0; ip < source->GetNbinsX(); ++ip)
  {
    for (int ir = 0; ir < source->GetNbinsY(); ++ir)
    {
      for (int iz = 0; iz < source->GetNbinsZ(); ++iz)
      {
        // do nothing if in TPOT acceptance
        bool in_range = mask->GetBinContent(ip + 1, ir + 1, iz + 1) > 0;
        if (in_range)
        {
          continue;
        }

        // phi not in TPOT range, rotate by steps of 2pi/12 (= one TPC sector) until found in range
        const double r = source->GetYaxis()->GetBinCenter(ir + 1);
        const double phi = source->GetXaxis()->GetBinCenter(ip + 1);
        static constexpr int n_sectors = 12;
        for (int sector = 1; sector < n_sectors; ++sector)
        {
          // get ref phi
          double phi_ref = phi + 2. * M_PI * sector / n_sectors;
          while (phi_ref >= 2 * M_PI)
          {
            phi_ref -= 2 * M_PI;
          };

          int ip_ref = source->GetXaxis()->FindBin(phi_ref) - 1;
          in_range = mask->GetBinContent(ip_ref + 1, ir + 1, iz + 1) > 0;
          if (!in_range)
          {
            continue;
          }

          // get normalization factor from CM histograms
          double scale = 1;
          if (source_cm)
          {
            const double distortion_local = source_cm->Interpolate(phi, r);
            const double distortion_ref = source_cm->Interpolate(phi_ref, r);
            scale = distortion_local / distortion_ref;
          }

          // update content and error
          const double distortion = scale * source->GetBinContent(ip_ref + 1, ir + 1, iz + 1);
          const double error = scale * source->GetBinError(ip_ref + 1, ir + 1, iz + 1);
          source->SetBinContent(ip + 1, ir + 1, iz + 1, distortion);
          source->SetBinError(ip + 1, ir + 1, iz + 1, error);
        }
      }
    }
  }
}

//_______________________________________________
void TpcSpaceChargeReconstructionHelper::extrapolate_phi2(TH3* source, const TH3* mask)
{
  if (!source)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_phi2 - invalid source histogram" << std::endl;
    return;
  }

  if (!mask)
  {
    std::cout << "TpcSpaceChargeReconstructionHelper::extrapolate_phi2 - invalid mask histogram" << std::endl;
    return;
  }

  // loop over phi bins
  for (int ir = 0; ir < source->GetNbinsY(); ++ir)
  {
    for (int iz = 0; iz < source->GetNbinsZ(); ++iz)
    {
      int ip_min = -1;
      for (int ip = 0; ip < source->GetNbinsX(); ++ip)
      {
        bool in_range = mask->GetBinContent(ip + 1, ir + 1, iz + 1) > 0;
        if (in_range)
        {
          ip_min = ip;
          continue;
        }

        // ip is not in range. Check if ip_min was set. If not, skip this bin. If yes, find next bin in range
        if (ip_min < 0)
        {
          continue;
        }

        int ip_max = ip + 1;
        for (; ip_max < source->GetNbinsX(); ++ip_max)
        {
          in_range = mask->GetBinContent(ip_max + 1, ir + 1, iz + 1) > 0;
          if (in_range)
          {
            break;
          }
        }

        // check that a valid bin was found
        if (ip_max == source->GetNbinsX())
        {
          continue;
        }

        // do the interpolation
        const double phi_min = source->GetXaxis()->GetBinCenter(ip_min + 1);
        const double phi_max = source->GetXaxis()->GetBinCenter(ip_max + 1);
        const double phi = source->GetXaxis()->GetBinCenter(ip + 1);

        const double alpha_min = (phi_max - phi) / (phi_max - phi_min);
        const double alpha_max = (phi - phi_min) / (phi_max - phi_min);

        const double distortion_min = source->GetBinContent(ip_min + 1, ir + 1, iz + 1);
        const double distortion_max = source->GetBinContent(ip_max + 1, ir + 1, iz + 1);
        const double distortion = alpha_min * distortion_min + alpha_max * distortion_max;

        const double error_min = source->GetBinError(ip_min + 1, ir + 1, iz + 1);
        const double error_max = source->GetBinError(ip_max + 1, ir + 1, iz + 1);
        const double error = std::sqrt(square(alpha_min * error_min) + square(alpha_max * error_max));

        source->SetBinContent(ip + 1, ir + 1, iz + 1, distortion);
        source->SetBinError(ip + 1, ir + 1, iz + 1, error);
      }
    }
  }
}

//_______________________________________________
std::tuple<TH3*, TH3*> TpcSpaceChargeReconstructionHelper::split(const TH3* source)
{
  if (!source)
  {
    return std::make_tuple<TH3*, TH3*>(nullptr, nullptr);
  }

  auto xaxis = source->GetXaxis();
  auto yaxis = source->GetYaxis();
  auto zaxis = source->GetZaxis();
  auto ibin = zaxis->FindBin((double) 0);

  // create histograms
  const TString name(source->GetName());
  auto hneg = new TH3F(
      name + "_negz", name + "_negz",
      xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax(),
      yaxis->GetNbins(), yaxis->GetXmin(), yaxis->GetXmax(),
      ibin - 1, zaxis->GetXmin(), zaxis->GetBinUpEdge(ibin - 1));

  auto hpos = new TH3F(
      name + "_posz", name + "_posz",
      xaxis->GetNbins(), xaxis->GetXmin(), xaxis->GetXmax(),
      yaxis->GetNbins(), yaxis->GetXmin(), yaxis->GetXmax(),
      zaxis->GetNbins() - (ibin - 1), zaxis->GetBinLowEdge(ibin), zaxis->GetXmax());

  // copy content and errors
  for (int ix = 0; ix < xaxis->GetNbins(); ++ix)
  {
    for (int iy = 0; iy < yaxis->GetNbins(); ++iy)
    {
      for (int iz = 0; iz < zaxis->GetNbins(); ++iz)
      {
        const auto content = source->GetBinContent(ix + 1, iy + 1, iz + 1);
        const auto error = source->GetBinError(ix + 1, iy + 1, iz + 1);

        if (iz < ibin - 1)
        {
          hneg->SetBinContent(ix + 1, iy + 1, iz + 1, content);
          hneg->SetBinError(ix + 1, iy + 1, iz + 1, error);
        }
        else
        {
          hpos->SetBinContent(ix + 1, iy + 1, iz - (ibin - 1) + 1, content);
          hpos->SetBinError(ix + 1, iy + 1, iz - (ibin - 1) + 1, error);
        }
      }
    }
  }

  // also copy axis titles
  for (const auto h : {hneg, hpos})
  {
    h->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
    h->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
    h->GetZaxis()->SetTitle(source->GetZaxis()->GetTitle());
  }

  return std::make_tuple(hneg, hpos);
}

//___________________________________________________________________________
TH3* TpcSpaceChargeReconstructionHelper::add_guarding_bins(const TH3* source, const TString& name)
{
  std::array<int, 3> bins{};
  std::array<double, 3> x_min{};
  std::array<double, 3> x_max{};

  int index = 0;
  for (const auto axis : {source->GetXaxis(), source->GetYaxis(), source->GetZaxis()})
  {
    // calculate bin width
    const auto bin_width = (axis->GetXmax() - axis->GetXmin()) / axis->GetNbins();

    // increase the number of bins by two
    bins[index] = axis->GetNbins() + 2;

    // update axis limits accordingly
    x_min[index] = axis->GetXmin() - bin_width;
    x_max[index] = axis->GetXmax() + bin_width;
    ++index;
  }

  // create new histogram
  auto hout = new TH3F(name, name,
                       bins[0], x_min[0], x_max[0],
                       bins[1], x_min[1], x_max[1],
                       bins[2], x_min[2], x_max[2]);

  // update axis legend
  hout->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
  hout->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
  hout->GetZaxis()->SetTitle(source->GetZaxis()->GetTitle());

  // copy content
  const auto phibins = source->GetXaxis()->GetNbins();
  const auto rbins = source->GetYaxis()->GetNbins();
  const auto zbins = source->GetZaxis()->GetNbins();

  // fill center
  for (int iphi = 0; iphi < phibins; ++iphi)
  {
    for (int ir = 0; ir < rbins; ++ir)
    {
      for (int iz = 0; iz < zbins; ++iz)
      {
        hout->SetBinContent(iphi + 2, ir + 2, iz + 2, source->GetBinContent(iphi + 1, ir + 1, iz + 1));
        hout->SetBinError(iphi + 2, ir + 2, iz + 2, source->GetBinError(iphi + 1, ir + 1, iz + 1));
      }
    }
  }

  // fill guarding phi bins
  /*
   * we use 2pi periodicity to do that:
   * - last valid bin is copied to first guarding bin;
   * - first valid bin is copied to last guarding bin
   */
  for (int ir = 0; ir < rbins + 2; ++ir)
  {
    for (int iz = 0; iz < zbins + 2; ++iz)
    {
      // copy last bin to first guarding bin
      hout->SetBinContent(1, ir + 1, iz + 1, hout->GetBinContent(phibins + 1, ir + 1, iz + 1));
      hout->SetBinError(1, ir + 1, iz + 1, hout->GetBinError(phibins + 1, ir + 1, iz + 1));

      // copy first bin to last guarding bin
      hout->SetBinContent(phibins + 2, ir + 1, iz + 1, hout->GetBinContent(2, ir + 1, iz + 1));
      hout->SetBinError(phibins + 2, ir + 1, iz + 1, hout->GetBinError(2, ir + 1, iz + 1));
    }
  }

  // fill guarding r bins
  for (int iphi = 0; iphi < phibins + 2; ++iphi)
  {
    for (int iz = 0; iz < zbins + 2; ++iz)
    {
      hout->SetBinContent(iphi + 1, 1, iz + 1, hout->GetBinContent(iphi + 1, 2, iz + 1));
      hout->SetBinError(iphi + 1, 1, iz + 1, hout->GetBinError(iphi + 1, 2, iz + 1));

      hout->SetBinContent(iphi + 1, rbins + 2, iz + 1, hout->GetBinContent(iphi + 1, rbins + 1, iz + 1));
      hout->SetBinError(iphi + 1, rbins + 2, iz + 1, hout->GetBinError(iphi + 1, rbins + 1, iz + 1));
    }
  }

  // fill guarding z bins
  for (int iphi = 0; iphi < phibins + 2; ++iphi)
  {
    for (int ir = 0; ir < rbins + 2; ++ir)
    {
      hout->SetBinContent(iphi + 1, ir + 1, 1, hout->GetBinContent(iphi + 1, ir + 1, 2));
      hout->SetBinError(iphi + 1, ir + 1, 1, hout->GetBinError(iphi + 1, ir + 1, 2));

      hout->SetBinContent(iphi + 1, ir + 1, zbins + 2, hout->GetBinContent(iphi + 1, ir + 1, zbins + 1));
      hout->SetBinError(iphi + 1, ir + 1, zbins + 2, hout->GetBinError(iphi + 1, ir + 1, zbins + 1));
    }
  }

  return hout;
}

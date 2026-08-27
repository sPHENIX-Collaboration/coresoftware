#include "PHGarfieldRossegger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TObject.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
  constexpr double epsilon0 = 8.8541878128e-12;
  constexpr double cmToM = 1.0e-2;
  constexpr double mToCm = 100.0;
  constexpr double elementaryCharge = 1.602176634e-19;
  constexpr unsigned int nTpcLayers = 48;
  constexpr unsigned int nTpcSectors = 12;
  constexpr unsigned int nTpcModules = 3;
  constexpr unsigned int nTpcSides = 2;
  constexpr unsigned int nLayersPerModule = 16;

  constexpr std::array<double, nTpcModules> firstLayerRadiusMm{{314.9835971194212, 416.5920253621078, 589.1096334932404}};
  constexpr std::array<double, nTpcModules> layerHeightMm{{5.657908935192669, 10.206891263554285, 10.970474549405283}};
  constexpr double r1GemInnerCm = 22.8;


  struct FrameDimensions
  {
    double inner_r_min_cm;
    double inner_r_max_cm;
    double outer_r_min_cm;
    double outer_r_max_cm;
    double left_width_cm;
    double right_width_cm;
  };

  constexpr unsigned int frameReferenceSector = 2;
  constexpr std::array<FrameDimensions, nTpcModules> frameDimensions{{
      {21.78, 22.20, 40.26, 40.61, 0.35, 0.52},
      {40.79, 41.14, 57.49, 57.84, 0.35, 0.52},
      {58.00, 58.35, 75.84, 76.28, 0.54, 0.54}}};
  constexpr double frameHalfAngle = 15.0 * std::numbers::pi / 180.0;



  unsigned int src_index(unsigned int ir, unsigned int ip, unsigned int iz, unsigned int np, unsigned int nz)
  { return (ir * np + ip) * nz + iz; }

  unsigned int fld_index(unsigned int ir, unsigned int ip, unsigned int iz, unsigned int np, unsigned int nz)
  { return (ir * np + ip) * nz + iz; }

  unsigned int mode_index(unsigned int in, unsigned int il, unsigned int nl)
  { return in * nl + il; }

  std::vector<double> edges(double lo, double hi, unsigned int n)
  {
    std::vector<double> out(n + 1);
    const double step = (hi - lo) / static_cast<double>(n);
    for (unsigned int i = 0; i <= n; ++i) { out[i] = lo + step * static_cast<double>(i); }
    return out;
  }

  std::vector<double> centers(const std::vector<double>& e)
  {
    std::vector<double> out(e.size() - 1);
    for (unsigned int i = 0; i + 1 < e.size(); ++i) { out[i] = 0.5 * (e[i] + e[i + 1]); }
    return out;
  }

  void append_edge_cm(std::vector<double>& edges_m, double value_cm)
  {
    edges_m.push_back(value_cm * cmToM);
  }

  void append_uniform_edges_cm(std::vector<double>& edges_m, double lo_cm, double hi_cm, double max_step_cm)
  {
    if (hi_cm <= lo_cm) { return; }
    if (edges_m.empty()) { append_edge_cm(edges_m, lo_cm); }
    const unsigned int n = std::max(1U, static_cast<unsigned int>(std::ceil((hi_cm - lo_cm) / max_step_cm)));
    for (unsigned int i = 1; i <= n; ++i)
    {
      append_edge_cm(edges_m, lo_cm + (hi_cm - lo_cm) * static_cast<double>(i) / static_cast<double>(n));
    }
  }

  void sort_unique_edges(std::vector<double>& edges_m)
  {
    std::sort(edges_m.begin(), edges_m.end());
    edges_m.erase(std::unique(edges_m.begin(), edges_m.end(), [](double lhs, double rhs) { return std::fabs(lhs - rhs) < 1.0e-12; }), edges_m.end());
  }

  std::vector<double> frame_observation_r_edges(double inner_cm, double outer_cm)
  {
    constexpr double fine_padding_cm = 1.20;
    std::vector<double> breakpoints{inner_cm, outer_cm};
    for (const FrameDimensions& f : frameDimensions)
    {
      for (const double boundary_cm : {f.inner_r_min_cm, f.inner_r_max_cm, f.outer_r_min_cm, f.outer_r_max_cm})
      {
        breakpoints.push_back(std::clamp(boundary_cm - fine_padding_cm, inner_cm, outer_cm));
        breakpoints.push_back(std::clamp(boundary_cm, inner_cm, outer_cm));
        breakpoints.push_back(std::clamp(boundary_cm + fine_padding_cm, inner_cm, outer_cm));
      }
    }
    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end(), [](double lhs, double rhs) { return std::fabs(lhs - rhs) < 1.0e-9; }), breakpoints.end());

    std::vector<double> out;
    for (unsigned int i = 0; i + 1 < breakpoints.size(); ++i)
    {
      const double lo_cm = breakpoints[i];
      const double hi_cm = breakpoints[i + 1];
      const double mid_cm = 0.5 * (lo_cm + hi_cm);
      bool fine = false;
      for (const FrameDimensions& f : frameDimensions)
      {
        fine = fine || (mid_cm >= f.inner_r_min_cm - fine_padding_cm && mid_cm <= f.inner_r_max_cm + fine_padding_cm) ||
                        (mid_cm >= f.outer_r_min_cm - fine_padding_cm && mid_cm <= f.outer_r_max_cm + fine_padding_cm);
      }
      append_uniform_edges_cm(out, lo_cm, hi_cm, fine ? 0.15 : 0.75);
    }
    sort_unique_edges(out);
    return out;
  }


  std::vector<double> frame_observation_z_edges(double half_length_cm)
  {
    const double near_padplane_cm = std::max(0.0, half_length_cm - 2.0);
    const double transition_cm = std::max(0.0, half_length_cm - 10.0);
    std::vector<double> out;
    append_uniform_edges_cm(out, 0.0, transition_cm, 5.0);
    append_uniform_edges_cm(out, transition_cm, near_padplane_cm, 0.50);
    append_uniform_edges_cm(out, near_padplane_cm, half_length_cm, 0.08);
    sort_unique_edges(out);
    return out;
  }


  double frame_sector_center(unsigned int side, unsigned int sector, double frame_reference_phi)
  {
    constexpr double sector_width = std::numbers::pi / 6.0;
    const int dsector = static_cast<int>(sector) - static_cast<int>(frameReferenceSector);
    double phi = side == 1U ? frame_reference_phi + dsector * sector_width : -frame_reference_phi - dsector * sector_width;
    const double twopi = 2.0 * std::numbers::pi;
    phi = std::fmod(phi, twopi);
    return phi < 0.0 ? phi + twopi : phi;
  }

  double local_frame_phi(unsigned int side, unsigned int sector, double global_phi, double frame_reference_phi)
  {
    const double center = frame_sector_center(side, sector, frame_reference_phi);
    const double delta = std::atan2(std::sin(global_phi - center), std::cos(global_phi - center));
    return side == 1U ? delta : -delta;
  }

  std::vector<double> frame_phi_edges(unsigned int, unsigned int side, double frame_reference_phi, double lo)
  {
    double fine_half_width = 0.0;
    for (const FrameDimensions& f : frameDimensions)
    {
      fine_half_width = std::max(fine_half_width, std::asin(f.left_width_cm / f.inner_r_max_cm));
      fine_half_width = std::max(fine_half_width, std::asin(f.right_width_cm / f.inner_r_max_cm));
    }
    constexpr double fine_step = 0.0015;
    constexpr double coarse_step = 0.05;
    const double hi = lo + 2.0 * std::numbers::pi;
    const double twopi = 2.0 * std::numbers::pi;

    auto shift_into_range = [lo, twopi](double phi)
    {
      while (phi < lo) { phi += twopi; }
      while (phi >= lo + twopi) { phi -= twopi; }
      return phi;
    };

    std::vector<double> breakpoints{lo, hi};
    auto add_edge = [&breakpoints, shift_into_range, lo, hi, fine_half_width](double phi)
    {
      const double folded = shift_into_range(phi);
      for (const double edge : {folded - fine_half_width, folded, folded + fine_half_width})
      {
        const double shifted = shift_into_range(edge);
        breakpoints.push_back(shifted);
        if (shifted < lo + fine_half_width) { breakpoints.push_back(shifted + 2.0 * std::numbers::pi); }
        if (shifted > hi - fine_half_width) { breakpoints.push_back(shifted - 2.0 * std::numbers::pi); }
      }
    };

    for (unsigned int sector = 0; sector < nTpcSectors; ++sector)
    {
      const double center = frame_sector_center(side, sector, frame_reference_phi);
      add_edge(center - frameHalfAngle);
      add_edge(center + frameHalfAngle);
    }

    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end(), [](double lhs, double rhs) { return std::fabs(lhs - rhs) < 1.0e-12; }), breakpoints.end());

    std::vector<double> out;
    for (unsigned int i = 0; i + 1 < breakpoints.size(); ++i)
    {
      const double plo = std::clamp(breakpoints[i], lo, hi);
      const double phi_hi = std::clamp(breakpoints[i + 1], lo, hi);
      if (phi_hi <= plo) { continue; }
      const double mid = 0.5 * (plo + phi_hi);
      bool fine = false;
      for (unsigned int sector = 0; sector < nTpcSectors; ++sector)
      {
        const double center = frame_sector_center(side, sector, frame_reference_phi);
        for (const double boundary : {center - frameHalfAngle, center + frameHalfAngle})
        {
          const double folded = shift_into_range(boundary);
          const double delta = std::atan2(std::sin(mid - folded), std::cos(mid - folded));
          fine = fine || std::fabs(delta) <= fine_half_width;
        }
      }
      if (out.empty()) { out.push_back(plo); }
      const unsigned int n = std::max(1U, static_cast<unsigned int>(std::ceil((phi_hi - plo) / (fine ? fine_step : coarse_step))));
      for (unsigned int j = 1; j <= n; ++j)
      {
        out.push_back(plo + (phi_hi - plo) * static_cast<double>(j) / static_cast<double>(n));
      }
    }
    sort_unique_edges(out);
    if (!out.empty()) { out.front() = lo; out.back() = hi; }
    return out;
  }



  std::vector<double> phi_edges(unsigned int n, double lo = 0.0)
  {
    std::vector<double> out(n + 1);
    const double step = 2.0 * std::numbers::pi / static_cast<double>(n);
    for (unsigned int i = 0; i <= n; ++i) { out[i] = lo + static_cast<double>(i) * step; }
    return out;
  }

  std::vector<double> phi_centers(unsigned int n, double lo = 0.0)
  {
    std::vector<double> out(n);
    const double step = 2.0 * std::numbers::pi / static_cast<double>(n);
    for (unsigned int i = 0; i < n; ++i) { out[i] = lo + (static_cast<double>(i) + 0.5) * step; }
    return out;
  }

  double wrap_phi(double phi)
  {
    const double twopi = 2.0 * std::numbers::pi;
    phi = std::fmod(phi, twopi);
    return phi < 0.0 ? phi + twopi : phi;
  }

  double wrapped_positive_delta(double phi_to, double phi_from)
  {
    return wrap_phi(phi_to - phi_from);
  }

  double wrapped_delta_phi(double phi_observation, double phi_source)
  {
    return std::atan2(std::sin(phi_observation - phi_source), std::cos(phi_observation - phi_source));
  }

  std::vector<std::pair<double, double>> periodic_segments(double start, double width)
  {
    const double twopi = 2.0 * std::numbers::pi;
    start = wrap_phi(start);
    const double end = start + width;
    if (end <= twopi + 1.0e-14) { return {{start, std::min(end, twopi)}}; }
    return {{start, twopi}, {0.0, end - twopi}};
  }

  double segment_overlap_length(const std::vector<std::pair<double, double>>& a, const std::vector<std::pair<double, double>>& b)
  {
    double total = 0.0;
    for (const auto& [alo, ahi] : a)
    {
      for (const auto& [blo, bhi] : b)
      {
        total += std::max(0.0, std::min(ahi, bhi) - std::max(alo, blo));
      }
    }
    return total;
  }

  double jprime(unsigned int m, double x)
  {
    return m == 0 ? -std::cyl_bessel_j(1, x) : 0.5 * (std::cyl_bessel_j(m - 1, x) - std::cyl_bessel_j(m + 1, x));
  }

  double yprime(unsigned int m, double x)
  {
    return m == 0 ? -std::cyl_neumann(1, x) : 0.5 * (std::cyl_neumann(m - 1, x) - std::cyl_neumann(m + 1, x));
  }

  double sinh_ratio(double k, double z, double len)
  {
    const double klen = k * len;
    if (klen < 50.0) { return std::sinh(k * z) / std::sinh(klen); }
    return std::exp(k * (z - len)) * (1.0 - std::exp(-2.0 * k * z)) / (1.0 - std::exp(-2.0 * klen));
  }

  double dz_sinh_ratio(double k, double z, double len)
  {
    const double klen = k * len;
    if (klen < 50.0) { return k * std::cosh(k * z) / std::sinh(klen); }
    return k * std::exp(k * (z - len)) * (1.0 + std::exp(-2.0 * k * z)) / (1.0 - std::exp(-2.0 * klen));
  }

  void scale_to_cm(std::vector<double>& values)
  {
    for (double& value : values) { value *= mToCm; }
  }
}  // namespace

PHGarfieldRossegger::PHGarfieldRossegger(const std::string& name)
  : SubsysReco(name)
{
}

void PHGarfieldRossegger::setGeometryCm(double inner_radius, double outer_radius, double half_length)
{
  m_aCm = inner_radius; m_bCm = outer_radius; m_lCm = half_length;
}

void PHGarfieldRossegger::setSourceRadiusCm(double min_radius, double max_radius)
{
  m_sourceRMinCm = min_radius; m_sourceRMaxCm = max_radius;
}

void PHGarfieldRossegger::setDensity(double reference_density_nC_per_m3, double k_eff, double radial_power_alpha)
{
  m_rhoReferenceNCPerM3 = reference_density_nC_per_m3; m_kEff = k_eff; m_radialPowerAlpha = radial_power_alpha;
}

void PHGarfieldRossegger::setPhiModulation(double m1_amplitude, double m1_phase, double m12_amplitude, double m12_phase)
{
  m_m1Amplitude = m1_amplitude; m_m1Phase = m1_phase; m_m12Amplitude = m12_amplitude; m_m12Phase = m12_phase;
}

void PHGarfieldRossegger::setSourceGrid(unsigned int nr, unsigned int nphi, unsigned int nz)
{
	std::cout<<"setSourceGrid nr nphi nz"<<nr<<" "<<nphi<<" "<<nz<<std::endl;
  m_nrSource = nr; m_nphiSource = nphi; m_nzSource = nz;
}

void PHGarfieldRossegger::setObservationGrid(unsigned int nr, unsigned int nphi, unsigned int nz)
{
  m_nrObs = nr; m_nphiObs = nphi; m_nzObs = nz;
}

void PHGarfieldRossegger::setModeTruncation(unsigned int m_phi_max, unsigned int n_radial_modes, unsigned int n_longitudinal_modes)
{
  m_mPhiMax = m_phi_max; m_nRadialModes = n_radial_modes; m_nLongitudinalModes = n_longitudinal_modes;
}

void PHGarfieldRossegger::setRadialJob(unsigned int job_index, unsigned int n_jobs)
{
  m_jobIndex = job_index; m_nJobs = std::max(1U, n_jobs);
}

void PHGarfieldRossegger::setGainHistograms(const std::string& side0_histogram, const std::string& side1_histogram)
{
  m_gainHistograms = {{side0_histogram, side1_histogram}};
}

void PHGarfieldRossegger::setDensityMapFile(const std::string& filename, const std::string& side0_histogram, const std::string& side1_histogram)
{
  m_densityMapFile = filename;
  m_densityMapHistograms = {{side0_histogram, side1_histogram}};
}

void PHGarfieldRossegger::setPHGarfieldField3DOutputFiles(const std::string& side0_filename, const std::string& side1_filename)
{
  m_phGarfieldField3DOutputFiles = {{side0_filename, side1_filename}};
}

int PHGarfieldRossegger::Init(PHCompositeNode*)
{
  return InitRun(nullptr);
}

int PHGarfieldRossegger::InitRun(PHCompositeNode*)
{
  if (m_done) { return Fun4AllReturnCodes::EVENT_OK; }
  const int status = calculate();
  m_done = status == Fun4AllReturnCodes::EVENT_OK;
  return status;
}

int PHGarfieldRossegger::process_event(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHGarfieldRossegger::validateConfig() const
{
  const double a = m_aCm * cmToM;
  const double b = m_bCm * cmToM;
  const double smin = m_sourceRMinCm * cmToM;
  const double smax = m_sourceRMaxCm * cmToM;
  const bool valid_source_radii = m_useDensityMap ?
      (0.0 < a && a <= smin && smin < smax && smax <= b) :
      (0.0 < a && a < smin && smin < smax && smax <= b);
  if (!m_useFrameChargeModel && !valid_source_radii)
  {
    std::cout << PHWHERE << " invalid Rossegger geometry/source radii" << std::endl;
    return false;
  }
  if (m_useDensityMap && (m_useFrameChargeModel || m_useRealTpcSourceGeometry || m_densityMapFile.empty()))
  {
    std::cout << PHWHERE << " density-map mode is incompatible with frame/real-TPC source modes or has no input file" << std::endl;
    return false;
  }
  if (m_lCm <= 0.0 || m_nrSource == 0 || m_nphiSource == 0 || m_nzSource == 0 ||
      m_nrObs == 0 || m_nphiObs == 0 || m_nzObs == 0 || m_nRadialModes == 0 || m_nLongitudinalModes == 0)
  {
    std::cout << PHWHERE << " invalid Rossegger grid or mode configuration" << std::endl;
    return false;
  }
  if (m_jobIndex >= m_nJobs)
  {
    std::cout << PHWHERE << " job index outside number of jobs" << std::endl;
    return false;
  }
  if (m_tpcSide >= nTpcSides)
  {
    std::cout << PHWHERE << " TPC side must be 0 or 1" << std::endl;
    return false;
  }
  return true;
}

unsigned int PHGarfieldRossegger::effectivePhiMax() const
{
  if (m_useDensityMap || m_useRealTpcSourceGeometry || m_useFrameChargeModel)
  {
    return m_mPhiMax;
  }
  const bool uniform = std::fabs(m_m1Amplitude) < 1.0e-15 && std::fabs(m_m12Amplitude) < 1.0e-15;
  return (m_autoAxisymmetric && uniform) ? 0 : m_mPhiMax;
}

std::pair<unsigned int, unsigned int> PHGarfieldRossegger::radialRange(unsigned int nr_obs) const
{
  return {(nr_obs * m_jobIndex) / m_nJobs, (nr_obs * (m_jobIndex + 1)) / m_nJobs};
}

double PHGarfieldRossegger::radialBoundaryFunction(double kval, unsigned int mode_m) const
{
  const double a = m_aCm * cmToM;
  const double b = m_bCm * cmToM;
  return std::cyl_bessel_j(mode_m, kval * b) * std::cyl_neumann(mode_m, kval * a) -
         std::cyl_neumann(mode_m, kval * b) * std::cyl_bessel_j(mode_m, kval * a);
}

double PHGarfieldRossegger::radialBasis(double kval, unsigned int mode_m, double radius_m) const
{
  const double a = m_aCm * cmToM;
  return std::cyl_bessel_j(mode_m, kval * radius_m) * std::cyl_neumann(mode_m, kval * a) -
         std::cyl_neumann(mode_m, kval * radius_m) * std::cyl_bessel_j(mode_m, kval * a);
}

double PHGarfieldRossegger::radialBasisDerivative(double kval, unsigned int mode_m, double radius_m) const
{
  const double a = m_aCm * cmToM;
  return kval * (jprime(mode_m, kval * radius_m) * std::cyl_neumann(mode_m, kval * a) -
                 yprime(mode_m, kval * radius_m) * std::cyl_bessel_j(mode_m, kval * a));
}

std::vector<double> PHGarfieldRossegger::findRadialRoots(unsigned int mode_m, unsigned int n_roots) const
{
  const double width = (m_bCm - m_aCm) * cmToM;
  const double step = 0.20 * std::numbers::pi / width;
  double left_k = std::max(1.0e-6, 0.10 * std::numbers::pi / width);
  double left_val = radialBoundaryFunction(left_k, mode_m);
  const double max_k = static_cast<double>(n_roots + mode_m + 20) * std::numbers::pi / width;
  std::vector<double> roots;

  for (double right_k = left_k + step; right_k <= max_k && roots.size() < n_roots; right_k += step)
  {
    const double right_val = radialBoundaryFunction(right_k, mode_m);
    if (std::isfinite(left_val) && std::isfinite(right_val) && left_val * right_val < 0.0)
    {
      double lo = right_k - step;
      double hi = right_k;
      double flo = radialBoundaryFunction(lo, mode_m);
      for (unsigned int iter = 0; iter < 120; ++iter)
      {
        const double mid = 0.5 * (lo + hi);
        const double fmid = radialBoundaryFunction(mid, mode_m);
        if (std::fabs(fmid) < 1.0e-14 || std::fabs(hi - lo) < 1.0e-12) { lo = mid; hi = mid; break; }
        if (flo * fmid <= 0.0) { hi = mid; }
        else { lo = mid; flo = fmid; }
      }
      const double root = 0.5 * (lo + hi);
      if (roots.empty() || std::fabs(root - roots.back()) > 1.0e-8) { roots.push_back(root); }
    }
    left_k = right_k;
    left_val = right_val;
  }
  if (roots.size() != n_roots)
  {
    throw std::runtime_error(std::format("Found {} roots for m={}, requested {}", roots.size(), mode_m, n_roots));
  }
  return roots;
}

std::vector<double> PHGarfieldRossegger::legendreRootsAndWeights(unsigned int n_points, std::vector<double>& weights) const
{
  std::vector<double> roots(n_points);
  weights.assign(n_points, 0.0);
  const unsigned int half = (n_points + 1) / 2;
  for (unsigned int i = 0; i < half; ++i)
  {
    double z = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n_points) + 0.5));
    double old_z = 0.0;
    double pp = 0.0;
    while (std::fabs(z - old_z) > 1.0e-15)
    {
      double p1 = 1.0;
      double p2 = 0.0;
      for (unsigned int j = 1; j <= n_points; ++j)
      {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * static_cast<double>(j) - 1.0) * z * p2 - (static_cast<double>(j) - 1.0) * p3) / static_cast<double>(j);
      }
      pp = static_cast<double>(n_points) * (z * p1 - p2) / (z * z - 1.0);
      old_z = z;
      z = old_z - p1 / pp;
    }
    roots[i] = -z;
    roots[n_points - 1 - i] = z;
    const double weight = 2.0 / ((1.0 - z * z) * pp * pp);
    weights[i] = weight;
    weights[n_points - 1 - i] = weight;
  }
  return roots;
}

PHGarfieldRossegger::SourceGrid PHGarfieldRossegger::buildSourceGrid() const
{
  if (m_useFrameChargeModel)
  {
    SourceGrid source_grid;
    std::vector<double>& redges = source_grid.r_edges_m;
    constexpr double max_step_cm = 0.25;

    auto add_interval = [&redges, max_step_cm](double rmin_cm, double rmax_cm)
    {
      redges.push_back(rmin_cm * cmToM);

      const double width_cm = rmax_cm - rmin_cm;
      const unsigned int n = std::max(1U, static_cast<unsigned int>(std::ceil(width_cm / max_step_cm)));
      for (unsigned int i = 1; i < n; ++i)
      {
        const double r_cm = rmin_cm + width_cm * static_cast<double>(i) / static_cast<double>(n);
        redges.push_back(r_cm * cmToM);
      }

      redges.push_back(rmax_cm * cmToM);
    };


    for (const FrameDimensions& f : frameDimensions)
    {
      add_interval(f.inner_r_min_cm, f.inner_r_max_cm);
      add_interval(f.inner_r_max_cm, f.outer_r_min_cm);
      add_interval(f.outer_r_min_cm, f.outer_r_max_cm);
    }
    std::sort(redges.begin(), redges.end());
    redges.erase(std::unique(redges.begin(), redges.end(), [](double lhs, double rhs) { return std::fabs(lhs - rhs) < 1.0e-12; }), redges.end());
    if (!(m_aCm * cmToM <= redges.front() && redges.back() <= m_bCm * cmToM))
    {
      throw std::runtime_error("Frame source radial grid is outside Rossegger field cage geometry");
    }
    source_grid.r_centers_m = centers(redges);
    source_grid.module_index.assign(source_grid.r_centers_m.size(), -1);
    source_grid.layer_index.assign(source_grid.r_centers_m.size(), -1);
    source_grid.is_antenna.assign(source_grid.r_centers_m.size(), false);
    for (unsigned int ir = 0; ir < source_grid.r_centers_m.size(); ++ir)
    {
      const double rcm = source_grid.r_centers_m[ir] * mToCm;
      for (unsigned int module = 0; module < frameDimensions.size(); ++module)
      {
        const FrameDimensions& f = frameDimensions[module];
        if (rcm >= f.inner_r_min_cm && rcm <= f.outer_r_max_cm)
        {
          source_grid.module_index[ir] = module;
          source_grid.layer_index[ir] = module;
          break;
        }
      }
    }
    return source_grid;
  }

  if (!m_useDensityMap && !m_useRealTpcSourceGeometry)
  {
    SourceGrid source_grid;
    source_grid.r_edges_m = edges(m_sourceRMinCm * cmToM, m_sourceRMaxCm * cmToM, m_nrSource);
    source_grid.r_centers_m = centers(source_grid.r_edges_m);
    source_grid.module_index.assign(source_grid.r_centers_m.size(), 0);
    source_grid.layer_index.assign(source_grid.r_centers_m.size(), 0);
    source_grid.is_antenna.assign(source_grid.r_centers_m.size(), false);
    return source_grid;
  }

  std::array<double, nTpcModules> first_cm{};
  std::array<double, nTpcModules> pitch_cm{};
  std::array<double, nTpcModules> first_low{};
  for (unsigned int module = 0; module < nTpcModules; ++module)
  {
    first_cm[module] = firstLayerRadiusMm[module] / 10.0;
    pitch_cm[module] = layerHeightMm[module] / 10.0;
    first_low[module] = first_cm[module] - 0.5 * pitch_cm[module];
  }

  SourceGrid source_grid;
  std::vector<double>& redges = source_grid.r_edges_m;

  // In real-TPC and density-map modes m_nrSource is intentionally unused; radial bins come from detector geometry.
  redges.push_back(r1GemInnerCm * cmToM);
  std::vector<double> inner_r1_edges_cm;
  for (double edge_cm = first_low[0] - pitch_cm[0]; edge_cm > r1GemInnerCm + 1.0e-12; edge_cm -= pitch_cm[0])
  {
    inner_r1_edges_cm.push_back(edge_cm);
  }
  for (auto iter = inner_r1_edges_cm.rbegin(); iter != inner_r1_edges_cm.rend(); ++iter)
  {
    redges.push_back(*iter * cmToM);
  }
  redges.push_back(first_low[0] * cmToM);

  for (unsigned int ilocal = 0; ilocal < nLayersPerModule; ++ilocal)
  {
    redges.push_back((first_cm[0] + (static_cast<double>(ilocal) + 0.5) * pitch_cm[0]) * cmToM);
  }

  for (unsigned int module = 1; module < nTpcModules; ++module)
  {
    redges.push_back(first_low[module] * cmToM);
    for (unsigned int ilocal = 0; ilocal < nLayersPerModule; ++ilocal)
    {
      redges.push_back((first_cm[module] + (static_cast<double>(ilocal) + 0.5) * pitch_cm[module]) * cmToM);
    }
  }

  for (unsigned int i = 1; i < redges.size(); ++i)
  {
    if (!(redges[i - 1] < redges[i])) { throw std::runtime_error("TPC source radial edges are not strictly increasing"); }
  }
  if (!(m_aCm * cmToM < redges.front() && redges.back() <= m_bCm * cmToM))
  {
    throw std::runtime_error("TPC source radial grid is outside Rossegger field cage geometry");
  }

  double min_width_cm = (redges[1] - redges[0]) * mToCm;
  double max_width_cm = min_width_cm;
  for (unsigned int i = 1; i < redges.size(); ++i)
  {
    const double width_cm = (redges[i] - redges[i - 1]) * mToCm;
    min_width_cm = std::min(min_width_cm, width_cm);
    max_width_cm = std::max(max_width_cm, width_cm);
  }
  std::cout << "  real-TPC source radial bins = " << (redges.size() - 1)
            << ", bin width min/max = " << min_width_cm << "/" << max_width_cm << " cm" << std::endl;

  source_grid.r_centers_m = centers(redges);
  const unsigned int nr = source_grid.r_centers_m.size();
  source_grid.module_index.assign(nr, -1);
  source_grid.layer_index.assign(nr, -1);
  source_grid.is_antenna.assign(nr, false);

  for (unsigned int ir = 0; ir < nr; ++ir)
  {
    const double rcm = source_grid.r_centers_m[ir] * mToCm;
    if (rcm >= r1GemInnerCm && rcm < first_low[0])
    {
      source_grid.module_index[ir] = 0;
      source_grid.layer_index[ir] = 0;
      source_grid.is_antenna[ir] = true;
      continue;
    }

    for (unsigned int module = 0; module < nTpcModules; ++module)
    {
      for (unsigned int ilocal = 0; ilocal < nLayersPerModule; ++ilocal)
      {
        const double lo = first_cm[module] + (static_cast<double>(ilocal) - 0.5) * pitch_cm[module];
        const double hi = first_cm[module] + (static_cast<double>(ilocal) + 0.5) * pitch_cm[module];
        if (rcm >= lo - 1.0e-12 && rcm < hi - 1.0e-12)
        {
          source_grid.module_index[ir] = module;
          source_grid.layer_index[ir] = module * nLayersPerModule + ilocal;
        }
      }
    }
  }

  return source_grid;
}

PHGarfieldRossegger::PadGeometry PHGarfieldRossegger::parsePadPlacementFile() const
{
  PadGeometry geometry{};
  std::array<std::array<std::array<bool, nTpcSectors>, nTpcSides>, nTpcModules> seen{};

  std::ifstream input(m_padPlacementFile);
  if (!input.is_open()) { throw std::runtime_error("Could not open pad placement file " + m_padPlacementFile); }

  const std::regex region_pattern(R"(region\s+(\d+)\s+first_layer\s+(\d+)\s+pads_per_sector\s+(\d+)\s+total_phibins\s+(\d+))");
  const std::regex entry_pattern(R"(side\s+([01])\s+sector\s+(\d+)\s+first_pad\s+(\d+)\s+first_phi\s+([-+0-9.eE]+)\s+last_pad\s+(\d+)\s+last_phi\s+([-+0-9.eE]+))");

  int current_region = -1;
  std::string line;
  std::smatch match;
  while (std::getline(input, line))
  {
    if (std::regex_search(line, match, region_pattern))
    {
      current_region = std::stoi(match[1].str());
      if (current_region < 0 || current_region >= static_cast<int>(nTpcModules))
      {
        throw std::runtime_error("Pad placement file has invalid region " + std::to_string(current_region));
      }
      continue;
    }

    if (current_region >= 0 && std::regex_search(line, match, entry_pattern))
    {
      const unsigned int side = std::stoul(match[1].str());
      const unsigned int sector = std::stoul(match[2].str());
      const int first_pad = std::stoi(match[3].str());
      const double first_phi = std::stod(match[4].str());
      const int last_pad = std::stoi(match[5].str());
      const double last_phi = std::stod(match[6].str());
      if (side >= nTpcSides || sector >= nTpcSectors) { throw std::runtime_error("Pad placement file has invalid side/sector"); }
      const int n_steps = last_pad - first_pad;
      if (n_steps <= 0) { throw std::runtime_error("Pad placement file has invalid pad indices"); }

      const double center_delta = wrapped_delta_phi(last_phi, first_phi);
      const double pad_pitch = center_delta / static_cast<double>(n_steps);
      const double start_edge = wrap_phi(first_phi - 0.5 * pad_pitch);
      const double end_edge = wrap_phi(last_phi + 0.5 * pad_pitch);
      const int direction = pad_pitch > 0.0 ? 1 : -1;
      const double angular_width = direction > 0 ? wrapped_positive_delta(end_edge, start_edge) : wrapped_positive_delta(start_edge, end_edge);

      geometry[current_region][side][sector] = {start_edge, angular_width, direction};
      seen[current_region][side][sector] = true;
    }
  }

  unsigned int nseen = 0;
  for (const auto& by_side : seen)
  {
    for (const auto& by_sector : by_side)
    {
      for (const bool value : by_sector) { if (value) { ++nseen; } }
    }
  }
  if (nseen != nTpcModules * nTpcSides * nTpcSectors)
  {
    throw std::runtime_error(std::format("Expected {} pad geometry entries, found {}", nTpcModules * nTpcSides * nTpcSectors, nseen));
  }

  return geometry;
}

PHGarfieldRossegger::GainMap PHGarfieldRossegger::readGainMap() const
{
  GainMap gain{};
  std::unique_ptr<TFile> input(TFile::Open(m_gainMapFile.c_str(), "READ"));
  if (!input || input->IsZombie()) { throw std::runtime_error("Could not open gain ROOT file " + m_gainMapFile); }

  for (unsigned int side = 0; side < nTpcSides; ++side)
  {
    auto* hist = dynamic_cast<TH2*>(input->Get(m_gainHistograms[side].c_str()));
    if (!hist) { throw std::runtime_error("Missing gain histogram " + m_gainHistograms[side] + " in " + m_gainMapFile); }
    if (hist->GetNbinsX() != static_cast<int>(nTpcSectors) || hist->GetNbinsY() != static_cast<int>(nTpcLayers))
    {
      throw std::runtime_error(std::format("{} has ({},{}) bins; expected ({},{})", m_gainHistograms[side], hist->GetNbinsX(), hist->GetNbinsY(), nTpcSectors, nTpcLayers));
    }
    for (unsigned int sector = 0; sector < nTpcSectors; ++sector)
    {
      for (unsigned int layer = 0; layer < nTpcLayers; ++layer)
      {
        const double value = hist->GetBinContent(sector + 1, layer + 1);
        if (!std::isfinite(value) || value <= 0.0)
        {
          throw std::runtime_error(std::format("Invalid gain value for side {}, sector {}, layer {}", side, sector, layer));
        }
        gain[side][sector][layer] = value;
      }
    }
  }
  return gain;
}

std::vector<double> PHGarfieldRossegger::makeChargeDensity(const SourceGrid& source_grid,
                                                           const std::vector<double>& phi_source_edges,
                                                           const std::vector<double>& z_source_edges_m,
                                                           unsigned int side) const
{
  const unsigned int nr = source_grid.r_centers_m.size();
  const unsigned int nphi = phi_source_edges.size() - 1;
  const unsigned int nz = z_source_edges_m.size() - 1;
  const double target_rho = m_kEff * m_rhoReferenceNCPerM3 * 1.0e-9;

  std::vector<double> volumes(nr * nphi * nz, 0.0);
  std::vector<double> raw(nr * nphi * nz, 0.0);
  double volume_sum = 0.0;

  if (m_useDensityMap)
  {
    std::unique_ptr<TFile> file(TFile::Open(m_densityMapFile.c_str(), "READ"));
    if (!file || file->IsZombie()) { throw std::runtime_error("Cannot open density-map file " + m_densityMapFile); }
    auto* hist = dynamic_cast<TH2*>(file->Get(m_densityMapHistograms[side].c_str()));
    if (!hist) { throw std::runtime_error("Missing density histogram " + m_densityMapHistograms[side]); }

    const auto ps = centers(phi_source_edges);
    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double period = xmax - xmin;

    for (unsigned int ir = 0; ir < nr; ++ir)
    {
      const double rvol = 0.5 * (source_grid.r_edges_m[ir + 1] * source_grid.r_edges_m[ir + 1] - source_grid.r_edges_m[ir] * source_grid.r_edges_m[ir]);
      const double rcm = source_grid.r_centers_m[ir] * mToCm;
      for (unsigned int ip = 0; ip < nphi; ++ip)
      {
        double phi = ps[ip];
        while (phi < xmin) { phi += period; }
        while (phi >= xmax) { phi -= period; }
        const double value = std::max(0.0, hist->Interpolate(phi, rcm));
        for (unsigned int iz = 0; iz < nz; ++iz)
        {
          const unsigned int idx = src_index(ir, ip, iz, nphi, nz);
          volumes[idx] = rvol * (phi_source_edges[ip + 1] - phi_source_edges[ip]) * (z_source_edges_m[iz + 1] - z_source_edges_m[iz]);
          raw[idx] = value;
          volume_sum += volumes[idx];
        }
      }
    }
  }
  else if (!m_useRealTpcSourceGeometry)
  {
    const auto ps = phi_centers(nphi);
    for (unsigned int ir = 0; ir < nr; ++ir)
    {
      const double rvol = 0.5 * (source_grid.r_edges_m[ir + 1] * source_grid.r_edges_m[ir + 1] - source_grid.r_edges_m[ir] * source_grid.r_edges_m[ir]);
      const double rshape = std::pow(source_grid.r_edges_m.front() / source_grid.r_centers_m[ir], m_radialPowerAlpha);
      for (unsigned int ip = 0; ip < nphi; ++ip)
      {
        const double pshape = 1.0 + m_m1Amplitude * std::cos(ps[ip] - m_m1Phase) + m_m12Amplitude * std::cos(12.0 * (ps[ip] - m_m12Phase));
        if (pshape < 0.0) { throw std::runtime_error("Phi modulation makes rho negative"); }
        for (unsigned int iz = 0; iz < nz; ++iz)
        {
          const unsigned int idx = src_index(ir, ip, iz, nphi, nz);
          volumes[idx] = rvol * (phi_source_edges[ip + 1] - phi_source_edges[ip]) * (z_source_edges_m[iz + 1] - z_source_edges_m[iz]);
          raw[idx] = rshape * pshape;
          volume_sum += volumes[idx];
        }
      }
    }
  }
  else
  {
    const PadGeometry pad_geometry = parsePadPlacementFile();
    const GainMap gain = readGainMap();
    const auto ps = phi_centers(nphi);
    std::vector<double> geometry_fraction(nr * nphi, 0.0);
    std::vector<double> inverse_gain_weight(nr * nphi, 0.0);

    for (unsigned int ir = 0; ir < nr; ++ir)
    {
      const int module = source_grid.module_index[ir];
      const int layer = source_grid.layer_index[ir];
      if (module < 0 || layer < 0) { continue; }

      for (unsigned int ip = 0; ip < nphi; ++ip)
      {
        const double bin_width = phi_source_edges[ip + 1] - phi_source_edges[ip];
        const auto bin_segments = periodic_segments(phi_source_edges[ip], bin_width);
        double fraction = 0.0;
        double weighted = 0.0;
        for (unsigned int sector = 0; sector < nTpcSectors; ++sector)
        {
          const PadInterval& interval = pad_geometry[module][side][sector];
          const double positive_start = interval.direction > 0 ? interval.start_edge : interval.start_edge - interval.angular_width;
          const auto active_segments = periodic_segments(positive_start, interval.angular_width);
          const double overlap = segment_overlap_length(bin_segments, active_segments) / bin_width;
          fraction += overlap;
          weighted += overlap * gain[side][sector][layer];
        }
        if (fraction > 1.0 + 1.0e-8)
        {
          throw std::runtime_error(std::format("Overlapping sector geometry for module {}, side {}, phi bin {}", module, side, ip));
        }
        geometry_fraction[ir * nphi + ip] = fraction;
        inverse_gain_weight[ir * nphi + ip] = weighted;
      }
    }

    double base_integral = 0.0;
    double corrected_integral = 0.0;
    for (unsigned int ir = 0; ir < nr; ++ir)
    {
      const double rvol = 0.5 * (source_grid.r_edges_m[ir + 1] * source_grid.r_edges_m[ir + 1] - source_grid.r_edges_m[ir] * source_grid.r_edges_m[ir]);
      const double rshape = std::pow(source_grid.r_edges_m.front() / source_grid.r_centers_m[ir], m_radialPowerAlpha);
      for (unsigned int ip = 0; ip < nphi; ++ip)
      {
        for (unsigned int iz = 0; iz < nz; ++iz)
        {
          const unsigned int idx = src_index(ir, ip, iz, nphi, nz);
          volumes[idx] = rvol * (phi_source_edges[ip + 1] - phi_source_edges[ip]) * (z_source_edges_m[iz + 1] - z_source_edges_m[iz]);
          //const double base_raw = rshape * geometry_fraction[ir * nphi + ip];
          //const double gain_raw = m_divideChargeByGain ? rshape * inverse_gain_weight[ir * nphi + ip] : base_raw;
	  const double pshape = 1.0 + m_m1Amplitude * std::cos(ps[ip] - m_m1Phase) + m_m12Amplitude * std::cos(12.0 * (ps[ip] - m_m12Phase));

	  if (pshape < 0.0)
	  {
	    throw std::runtime_error("Phi modulation makes rho negative");
	  }

	  const double base_raw = rshape * pshape * geometry_fraction[ir * nphi + ip];
	  const double gain_raw = m_divideChargeByGain ? rshape * pshape * inverse_gain_weight[ir * nphi + ip] : base_raw;

	  raw[idx] = gain_raw;
          base_integral += base_raw * volumes[idx];
          corrected_integral += gain_raw * volumes[idx];
          volume_sum += volumes[idx];
        }
      }
    }

    if (m_divideChargeByGain && m_normalizeGainWeightedTotal)
    {
      if (base_integral <= 0.0 || corrected_integral <= 0.0) { throw std::runtime_error("Invalid gain-weighted source normalization"); }
      const double gain_normalization = corrected_integral / base_integral;
      for (double& value : raw) { value /= gain_normalization; }
      std::cout << "  gain normalization side " << side << " = " << gain_normalization << std::endl;
    }
  }

  double shape_sum = 0.0;
  for (unsigned int idx = 0; idx < raw.size(); ++idx) { shape_sum += raw[idx] * volumes[idx]; }
  const double shape_mean = shape_sum / volume_sum;
  if (shape_mean <= 0.0) { throw std::runtime_error("Charge-density shape has zero volume average"); }

  const double scale = (m_useDensityMap && !m_normalizeDensityMap) ? 1.0e-9 : target_rho / shape_mean;
  std::vector<double> rho(raw.size(), 0.0);
  double total_charge = 0.0;
  for (unsigned int idx = 0; idx < rho.size(); ++idx)
  {
    rho[idx] = raw[idx] * scale;
    total_charge += rho[idx] * volumes[idx];
  }

  if (m_useDensityMap)
  {
    std::cout << "  density map side " << side << ": " << m_densityMapHistograms[side]
              << ", input mean = " << shape_mean
              << (m_normalizeDensityMap ? " arb., scale = " : " nC/m^3, scale = ") << scale << std::endl;
  }
  std::cout << "  mean rho side " << side << " = " << total_charge / volume_sum
            << " C/m^3, charge = " << total_charge << " C, ions = " << total_charge / elementaryCharge << std::endl;
  return rho;
}


bool PHGarfieldRossegger::pointInPolygon(const std::vector<Point2D>& polygon, double x_cm, double y_cm)
{
  bool inside = false;
  const std::size_t n = polygon.size();
  if (n < 3) { return false; }

  for (std::size_t i = 0, j = n - 1; i < n; j = i++)
  {
    const double xi = polygon[i].x_cm;
    const double yi = polygon[i].y_cm;
    const double xj = polygon[j].x_cm;
    const double yj = polygon[j].y_cm;
    const bool intersects = ((yi > y_cm) != (yj > y_cm)) &&
                            (x_cm < (xj - xi) * (y_cm - yi) / (yj - yi) + xi);
    if (intersects) { inside = !inside; }
  }
  return inside;
}

double PHGarfieldRossegger::polygonAreaCm2(const std::vector<Point2D>& polygon)
{
  const std::size_t n = polygon.size();
  if (n < 3) { return 0.0; }

  double area = 0.0;
  for (std::size_t i = 0, j = n - 1; i < n; j = i++)
  {
    area += polygon[j].x_cm * polygon[i].y_cm - polygon[i].x_cm * polygon[j].y_cm;
  }
  return 0.5 * std::fabs(area);
}

void PHGarfieldRossegger::loadFramePolygons() const
{
  if (m_framePolygonsLoaded) { return; }

  FramePolygons polygons{};
  std::ifstream input(m_frameGeometryFile);
  if (!input.is_open()) { throw std::runtime_error("Could not open frame geometry CSV " + m_frameGeometryFile); }

  std::string line;
  std::getline(input, line);
  while (std::getline(input, line))
  {
    if (line.empty()) { continue; }

    std::stringstream ss(line);
    std::array<std::string, 5> fields{};
    for (std::string& field : fields)
    {
      if (!std::getline(ss, field, ',')) { throw std::runtime_error("Malformed frame geometry CSV row: " + line); }
    }

    const unsigned int module = std::stoul(fields[0]);
    if (module >= nTpcModules) { throw std::runtime_error("Frame geometry CSV has invalid module " + fields[0]); }

    std::vector<Point2D>* polygon = nullptr;
    if (fields[1] == "inner") { polygon = &polygons[module].inner; }
    else if (fields[1] == "outer") { polygon = &polygons[module].outer; }
    else { throw std::runtime_error("Frame geometry CSV has invalid boundary " + fields[1]); }

    polygon->push_back({std::stod(fields[3]) / 10.0, std::stod(fields[4]) / 10.0});
  }

  for (unsigned int module = 0; module < nTpcModules; ++module)
  {
    if (polygons[module].inner.size() < 3 || polygons[module].outer.size() < 3)
    {
      throw std::runtime_error(std::format("Frame geometry CSV missing inner/outer polygon for module {}", module));
    }
  }

  m_framePolygons = std::move(polygons);
  m_framePolygonsLoaded = true;
}

bool PHGarfieldRossegger::pointInFrameGeometry(const FramePolygon& frame, double r_cm, double phi_rel) const
{
  const double x_cm = r_cm * std::sin(phi_rel);
  const double y_cm = r_cm * std::cos(phi_rel);
  return pointInPolygon(frame.outer, x_cm, y_cm) && !pointInPolygon(frame.inner, x_cm, y_cm);
}

double PHGarfieldRossegger::frameAreaM2(const FramePolygon& frame) const
{
  return (polygonAreaCm2(frame.outer) - polygonAreaCm2(frame.inner)) * 1.0e-4;
}

PHGarfieldRossegger::FrameBoundaryPattern PHGarfieldRossegger::makeFrameBoundaryPattern(const SourceGrid& source_grid,
                                                                                         const std::vector<double>& phi_source_edges,
                                                                                         unsigned int side) const
{
  const unsigned int nr = source_grid.r_centers_m.size();
  const unsigned int nphi = phi_source_edges.size() - 1;
  FrameBoundaryPattern pattern;
  pattern.geometry_fraction.assign(nr * nphi, 0.0);
  pattern.weight.assign(nr * nphi, 0.0);
  pattern.boundary_potential.assign(nr * nphi, 0.0);

  loadFramePolygons();
  std::array<double, nTpcModules> frame_area_by_module_m2{};
  double reference_piece_area = 0.0;
  for (unsigned int module = 0; module < nTpcModules; ++module)
  {
    const double inner_area_m2 = polygonAreaCm2(m_framePolygons[module].inner) * 1.0e-4;
    const double outer_area_m2 = polygonAreaCm2(m_framePolygons[module].outer) * 1.0e-4;
    frame_area_by_module_m2[module] = outer_area_m2 - inner_area_m2;
    if (frame_area_by_module_m2[module] <= 0.0)
    {
      throw std::runtime_error(std::format("Frame geometry CSV gives non-positive area for module {}", module));
    }
    std::cout << "  frame geometry module " << module
              << ": outer area = " << outer_area_m2
              << " m^2, inner opening area = " << inner_area_m2
              << " m^2, frame area = " << frame_area_by_module_m2[module]
              << " m^2" << std::endl;
    reference_piece_area += frame_area_by_module_m2[module];
  }
  reference_piece_area /= static_cast<double>(nTpcModules);

  constexpr unsigned int n_r_samples = 5;
  constexpr unsigned int n_phi_samples = 5;
  constexpr double inv_samples = 1.0 / static_cast<double>(n_r_samples * n_phi_samples);
  double weighted_area = 0.0;
  for (unsigned int ir = 0; ir < nr; ++ir)
  {
    const int module = source_grid.module_index[ir];
    if (module < 0) { continue; }
    const double r1_m = source_grid.r_edges_m[ir];
    const double r2_m = source_grid.r_edges_m[ir + 1];
    const double r1_cm = r1_m * mToCm;
    const double r2_cm = r2_m * mToCm;
    const double r2_span_cm2 = r2_cm * r2_cm - r1_cm * r1_cm;
    const double area_r = 0.5 * (r2_m * r2_m - r1_m * r1_m);

    for (unsigned int ip = 0; ip < nphi; ++ip)
    {
      const double bin_width = phi_source_edges[ip + 1] - phi_source_edges[ip];
      double geometry_fraction = 0.0;
      double weight = 0.0;

      for (unsigned int irs = 0; irs < n_r_samples; ++irs)
      {
        const double radial_fraction = (static_cast<double>(irs) + 0.5) / static_cast<double>(n_r_samples);
        const double r_cm = std::sqrt(r1_cm * r1_cm + radial_fraction * r2_span_cm2);
        for (unsigned int ips = 0; ips < n_phi_samples; ++ips)
        {
          const double phi = phi_source_edges[ip] + (static_cast<double>(ips) + 0.5) * bin_width / static_cast<double>(n_phi_samples);
          for (unsigned int sector = 0; sector < nTpcSectors; ++sector)
          {
            const double phi_rel = local_frame_phi(side, sector, phi, m_frameReferencePhi);
            if (!pointInFrameGeometry(m_framePolygons[static_cast<unsigned int>(module)], r_cm, phi_rel)) { continue; }
            geometry_fraction += inv_samples;
            if (m_frameChargeWeighting == FrameChargeWeighting::EqualChargePerPiece)
            {
              weight += inv_samples * reference_piece_area / frame_area_by_module_m2[module];
            }
            else
            {
              weight += inv_samples;
            }
          }
        }
      }

      const unsigned int idx = ir * nphi + ip;
      pattern.geometry_fraction[idx] = std::clamp(geometry_fraction, 0.0, 1.0);
      pattern.weight[idx] = weight;
      pattern.boundary_potential[idx] = m_frameBoundaryPotential * weight;
      weighted_area += area_r * bin_width * weight;
    }
  }

  std::cout << "  frame boundary potential side " << side << " = " << m_frameBoundaryPotential
            << " V, weighting = "
            << (m_frameChargeWeighting == FrameChargeWeighting::EqualChargePerPiece ? "EqualChargePerPiece" : "ProportionalToArea")
            << ", frame geometry = " << m_frameGeometryFile
            << ", frame reference sector = " << frameReferenceSector
            << ", weighted transverse area = " << weighted_area << " m^2" << std::endl;
  return pattern;
}

int PHGarfieldRossegger::calculate()
{
  if (!validateConfig()) { return Fun4AllReturnCodes::ABORTRUN; }
  try
  {
    const double a = m_aCm * cmToM;
    const double b = m_bCm * cmToM;
    const double len = m_lCm * cmToM;
    const unsigned int mmax = effectivePhiMax();

    std::cout << Name() << " Rossegger field calculation" << std::endl;
    if (m_useFrameChargeModel)
    {
      std::cout << "  using frame boundary-potential source geometry" << std::endl;
    }
    else if (m_useDensityMap)
    {
      std::cout << "  using ROOT r-phi density map: " << m_densityMapFile
                << (m_normalizeDensityMap ? " (normalized to setDensity reference)" : " (input values in nC/m^3)") << std::endl;
    }
    else if (m_useRealTpcSourceGeometry)
    {
      std::cout << "  using real TPC radial source geometry and layer/sector gain maps" << std::endl;
    }

    const SourceGrid source_grid = buildSourceGrid();
    const unsigned int nr_source = source_grid.r_centers_m.size();
    auto pse = m_useFrameChargeModel ? frame_phi_edges(m_nphiSource, m_tpcSide, m_frameReferencePhi, -std::numbers::pi) : ((m_useDensityMap || m_useRealTpcSourceGeometry) ? phi_edges(m_nphiSource, -std::numbers::pi) : phi_edges(m_nphiSource));
    auto zse = edges(0.0, len, m_nzSource);
    auto roe = m_useFrameChargeModel ? frame_observation_r_edges(m_aCm, m_bCm) : edges(a, b, m_nrObs);
    auto poe = m_useFrameChargeModel ? frame_phi_edges(m_nphiObs, m_tpcSide, m_frameReferencePhi, -std::numbers::pi) : ((m_useDensityMap || m_writePHGarfieldField3D) ? phi_edges(m_nphiObs, -std::numbers::pi) : phi_edges(m_nphiObs));
    auto zoe = m_useFrameChargeModel ? frame_observation_z_edges(m_lCm) : edges(0.0, len, m_nzObs);
    m_nphiSource = static_cast<unsigned int>(pse.size() - 1);
    m_nphiObs = static_cast<unsigned int>(poe.size() - 1);
    m_nzSource = static_cast<unsigned int>(zse.size() - 1);
    m_nrObs = static_cast<unsigned int>(roe.size() - 1);
    m_nzObs = static_cast<unsigned int>(zoe.size() - 1);
    const auto [r_begin, r_end] = radialRange(m_nrObs);
    std::cout << "  radial job " << m_jobIndex << "/" << m_nJobs << " bins [" << r_begin << ", " << r_end << ")" << std::endl;
    const auto ps = centers(pse);
    const auto zs = centers(zse);
    const auto ro = centers(roe);
    const auto po = centers(poe);
    const auto zo = centers(zoe);

    std::vector<std::vector<double>> km(mmax + 1);
    std::vector<std::vector<double>> norms(mmax + 1);
    std::vector<double> wg;
    const auto xg = legendreRootsAndWeights(192, wg);
    std::vector<double> rg(xg.size(), 0.0);
    const double jac = 0.5 * (b - a);
    for (unsigned int ig = 0; ig < xg.size(); ++ig) { rg[ig] = jac * xg[ig] + 0.5 * (b + a); }
    for (unsigned int im = 0; im <= mmax; ++im)
    {
      km[im] = findRadialRoots(im, m_nRadialModes);
      norms[im].assign(m_nRadialModes, 0.0);
      for (unsigned int in = 0; in < m_nRadialModes; ++in)
      {
        double sum = 0.0;
        for (unsigned int ig = 0; ig < xg.size(); ++ig)
        {
          const double rb = radialBasis(km[im][in], im, rg[ig]);
          sum += wg[ig] * rg[ig] * rb * rb;
        }
        norms[im][in] = jac * sum;
      }
    }

    std::array<bool, nTpcSides> run_side{{false, false}};
    run_side[m_tpcSide] = true;
    if (m_writePHGarfieldField3D)
    {
      run_side[0] = true;
      run_side[1] = true;
    }

    if (m_useFrameChargeModel)
    {
      std::vector<double> wg_source;
      const auto xg_source = legendreRootsAndWeights(8, wg_source);
      for (unsigned int side = 0; side < nTpcSides; ++side)
      {
        if (!run_side[side]) { continue; }
        const FrameBoundaryPattern frame_pattern = makeFrameBoundaryPattern(source_grid, pse, side);

        std::vector<std::vector<double>> mc(mmax + 1), ms(mmax + 1);
        for (unsigned int im = 0; im <= mmax; ++im)
        {
          const double anorm = im == 0 ? 2.0 * std::numbers::pi : std::numbers::pi;
          mc[im].assign(m_nRadialModes, 0.0);
          ms[im].assign(m_nRadialModes, 0.0);
          for (unsigned int in = 0; in < m_nRadialModes; ++in)
          {
            double cproj = 0.0;
            double sproj = 0.0;
            for (unsigned int irs = 0; irs < nr_source; ++irs)
            {
              const double r1 = source_grid.r_edges_m[irs];
              const double r2 = source_grid.r_edges_m[irs + 1];
              const double radial_jac = 0.5 * (r2 - r1);
              const double radial_center = 0.5 * (r2 + r1);
              double radial_integral = 0.0;
              for (unsigned int ig = 0; ig < xg_source.size(); ++ig)
              {
                const double r = radial_jac * xg_source[ig] + radial_center;
                const double rb = radialBasis(km[im][in], im, r);
                radial_integral += wg_source[ig] * r * rb;
              }
              radial_integral *= radial_jac;
              for (unsigned int ips = 0; ips < m_nphiSource; ++ips)
              {
                const unsigned int idx = irs * m_nphiSource + ips;
                const double boundary_potential = frame_pattern.boundary_potential[idx];
                if (boundary_potential == 0.0) { continue; }
                const double phi1 = pse[ips];
                const double phi2 = pse[ips + 1];
                double int_cos = 0.0;
                double int_sin = 0.0;
                if (im == 0)
                {
                  int_cos = phi2 - phi1;
                }
                else
                {
                  const double m = static_cast<double>(im);
                  int_cos = (std::sin(m * phi2) - std::sin(m * phi1)) / m;
                  int_sin = (std::cos(m * phi1) - std::cos(m * phi2)) / m;
                }
                cproj += boundary_potential * int_cos * radial_integral;
                if (im > 0) { sproj += boundary_potential * int_sin * radial_integral; }
              }
            }
            mc[im][in] = cproj / (norms[im][in] * anorm);
            ms[im][in] = im == 0 ? 0.0 : sproj / (norms[im][in] * anorm);
          }
        }
        for (unsigned int in = 0; in < std::min(3U, m_nRadialModes); ++in)
        {
          mc[0][in] = 0.0;
        }

        const unsigned int fsize = m_nrObs * m_nphiObs * m_nzObs;
        std::vector<double> phi(fsize, 0.0), er(fsize, 0.0), ep(fsize, 0.0), ez(fsize, 0.0);
        for (unsigned int iro = r_begin; iro < r_end; ++iro)
        {
          for (unsigned int im = 0; im <= mmax; ++im)
          {
            std::vector<double> rb(m_nRadialModes, 0.0), drb(m_nRadialModes, 0.0);
            for (unsigned int in = 0; in < m_nRadialModes; ++in)
            {
              rb[in] = radialBasis(km[im][in], im, ro[iro]);
              drb[in] = radialBasisDerivative(km[im][in], im, ro[iro]);
            }
            for (unsigned int ipo = 0; ipo < m_nphiObs; ++ipo)
            {
              const double cp = std::cos(static_cast<double>(im) * po[ipo]);
              const double sp = std::sin(static_cast<double>(im) * po[ipo]);
              for (unsigned int izo = 0; izo < m_nzObs; ++izo)
              {
                const unsigned int idx = fld_index(iro, ipo, izo, m_nphiObs, m_nzObs);
                for (unsigned int in = 0; in < m_nRadialModes; ++in)
                {
                  const double kval = km[im][in];
                  const double ac = mc[im][in];
                  const double as = ms[im][in];
                  const double amp = ac * cp + as * sp;
                  const double zfac = sinh_ratio(kval, zo[izo], len);
                  const double dzfac = dz_sinh_ratio(kval, zo[izo], len);
                  phi[idx] += rb[in] * zfac * amp;
                  er[idx] -= drb[in] * zfac * amp;
                  ez[idx] -= m_frameEzScale * rb[in] * dzfac * amp;
                  if (im > 0) { ep[idx] += (static_cast<double>(im) / ro[iro]) * rb[in] * zfac * (ac * sp - as * cp); }
                }
              }
            }
          }
        }

        if (side == m_tpcSide)
        {
          writeGarfieldRootFile(roe, poe, zoe, er, ep, ez, r_begin, r_end);
          if (m_writeField3D) { writeFrameBoundaryField3DRootFile(source_grid.r_edges_m, roe, pse, poe, zoe, frame_pattern, phi, er, ep, ez, r_begin, r_end); }
        }
        if (m_writePHGarfieldField3D) { writePHGarfieldField3DRootFile(roe, poe, zoe, phi, er, ep, ez, side, r_begin, r_end); }
      }
    }
    else
    {
      std::vector<double> q(m_nLongitudinalModes, 0.0);
      for (unsigned int il = 0; il < m_nLongitudinalModes; ++il) { q[il] = static_cast<double>(il + 1) * std::numbers::pi / len; }

      for (unsigned int side = 0; side < nTpcSides; ++side)
      {
        if (!run_side[side]) { continue; }
        const std::vector<double> rho = makeChargeDensity(source_grid, pse, zse, side);

        std::vector<std::vector<double>> mc(mmax + 1), ms(mmax + 1);
        for (unsigned int im = 0; im <= mmax; ++im)
        {
          mc[im].assign(m_nRadialModes * m_nLongitudinalModes, 0.0);
          ms[im].assign(m_nRadialModes * m_nLongitudinalModes, 0.0);
          for (unsigned int in = 0; in < m_nRadialModes; ++in)
          {
            for (unsigned int il = 0; il < m_nLongitudinalModes; ++il)
            {
              double cproj = 0.0;
              double sproj = 0.0;
              for (unsigned int irs = 0; irs < nr_source; ++irs)
              {
                const double rb = radialBasis(km[im][in], im, source_grid.r_centers_m[irs]);
                const double rvol = 0.5 * (source_grid.r_edges_m[irs + 1] * source_grid.r_edges_m[irs + 1] - source_grid.r_edges_m[irs] * source_grid.r_edges_m[irs]);
                for (unsigned int ips = 0; ips < m_nphiSource; ++ips)
                {
                  const double cp = std::cos(static_cast<double>(im) * ps[ips]);
                  const double sp = std::sin(static_cast<double>(im) * ps[ips]);
                  const double dphi = pse[ips + 1] - pse[ips];
                  for (unsigned int izs = 0; izs < m_nzSource; ++izs)
                  {
                    const unsigned int idx = src_index(irs, ips, izs, m_nphiSource, m_nzSource);
                    const double dz = zse[izs + 1] - zse[izs];
                    const double charge = rho[idx] * rvol * dphi * dz;
                    const double zfactor = std::sin(zs[izs] * q[il]);
                    cproj += charge * rb * cp * zfactor;
                    sproj += charge * rb * sp * zfactor;
                  }
                }
              }
              mc[im][mode_index(in, il, m_nLongitudinalModes)] = cproj;
              ms[im][mode_index(in, il, m_nLongitudinalModes)] = im == 0 ? 0.0 : sproj;
            }
          }
        }

        const unsigned int fsize = m_nrObs * m_nphiObs * m_nzObs;
        std::vector<double> phi(fsize, 0.0), er(fsize, 0.0), ep(fsize, 0.0), ez(fsize, 0.0);
        for (unsigned int iro = r_begin; iro < r_end; ++iro)
        {
          for (unsigned int im = 0; im <= mmax; ++im)
          {
            const double anorm = im == 0 ? 2.0 * std::numbers::pi : std::numbers::pi;
            std::vector<double> rb(m_nRadialModes, 0.0), drb(m_nRadialModes, 0.0);
            for (unsigned int in = 0; in < m_nRadialModes; ++in)
            {
              rb[in] = radialBasis(km[im][in], im, ro[iro]);
              drb[in] = radialBasisDerivative(km[im][in], im, ro[iro]);
            }
            for (unsigned int ipo = 0; ipo < m_nphiObs; ++ipo)
            {
              const double cp = std::cos(static_cast<double>(im) * po[ipo]);
              const double sp = std::sin(static_cast<double>(im) * po[ipo]);
              for (unsigned int izo = 0; izo < m_nzObs; ++izo)
              {
                const unsigned int idx = fld_index(iro, ipo, izo, m_nphiObs, m_nzObs);
                for (unsigned int in = 0; in < m_nRadialModes; ++in)
                {
                  for (unsigned int il = 0; il < m_nLongitudinalModes; ++il)
                  {
                    const double kval = km[im][in];
                    const double qval = q[il];
                    const unsigned int midx = mode_index(in, il, m_nLongitudinalModes);
                    const double denom = epsilon0 * norms[im][in] * (len / 2.0) * anorm * (kval * kval + qval * qval);
                    const double ac = mc[im][midx] / denom;
                    const double as = ms[im][midx] / denom;
                    const double amp = ac * cp + as * sp;
                    const double sz = std::sin(zo[izo] * qval);
                    const double qcz = qval * std::cos(zo[izo] * qval);
                    phi[idx] += rb[in] * sz * amp;
                    er[idx] -= drb[in] * sz * amp;
                    ez[idx] -= rb[in] * qcz * amp;
                    if (im > 0) { ep[idx] += (static_cast<double>(im) / ro[iro]) * rb[in] * sz * (ac * sp - as * cp); }
                  }
                }
              }
            }
          }
        }

        if (side == m_tpcSide)
        {
          writeGarfieldRootFile(roe, poe, zoe, er, ep, ez, r_begin, r_end);
          if (m_writeField3D) { writeField3DRootFile(source_grid.r_edges_m, roe, pse, poe, zse, zoe, rho, phi, er, ep, ez, r_begin, r_end); }
        }
        if (m_writePHGarfieldField3D) { writePHGarfieldField3DRootFile(roe, poe, zoe, phi, er, ep, ez, side, r_begin, r_end); }
      }
    }

    if (m_verifyOutput && !verifyOutput()) { return Fun4AllReturnCodes::ABORTRUN; }
  }
  catch (const std::exception& e)
  {
    std::cout << PHWHERE << " " << Name() << " failed: " << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfieldRossegger::writeGarfieldRootFile(const std::vector<double>& r_edges_m,
                                                const std::vector<double>& phi_edges,
                                                const std::vector<double>& z_edges_m,
                                                const std::vector<double>& er,
                                                const std::vector<double>& ephi,
                                                const std::vector<double>& ez,
                                                unsigned int r_begin,
                                                unsigned int r_end) const
{
  auto rcm = r_edges_m;
  auto zcm = z_edges_m;
  scale_to_cm(rcm);
  scale_to_cm(zcm);

  TFile file(m_garfieldOutputFile.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + m_garfieldOutputFile); }
  file.mkdir("QA");
  file.cd("QA");
  TH2D hEr("hErDefault", "Axisymmetric radial field for PHGarfield", m_nrObs, rcm.data(), m_nzObs, zcm.data());
  TH2D hEphi("hEphiDefault", "Axisymmetric azimuthal field for PHGarfield", m_nrObs, rcm.data(), m_nzObs, zcm.data());
  TH2D hEz("hEzDefault", "Axisymmetric local longitudinal field for PHGarfield", m_nrObs, rcm.data(), m_nzObs, zcm.data());
  TH2D hErRPhi("hErRPhi", "Radial field for PHGarfield averaged over |z|;r [cm];#phi [rad]", m_nrObs, rcm.data(), m_nphiObs, phi_edges.data());
  TH2D hEphiRPhi("hEphiRPhi", "Azimuthal field for PHGarfield averaged over |z|;r [cm];#phi [rad]", m_nrObs, rcm.data(), m_nphiObs, phi_edges.data());
  TH2D hEzRPhi("hEzRPhi", "Local longitudinal field for PHGarfield averaged over |z|;r [cm];#phi [rad]", m_nrObs, rcm.data(), m_nphiObs, phi_edges.data());
  TH2D hErZPhi("hErZPhi", "Radial field for PHGarfield averaged over r;|z| [cm];#phi [rad]", m_nzObs, zcm.data(), m_nphiObs, phi_edges.data());
  TH2D hEphiZPhi("hEphiZPhi", "Azimuthal field for PHGarfield averaged over r;|z| [cm];#phi [rad]", m_nzObs, zcm.data(), m_nphiObs, phi_edges.data());
  TH2D hEzZPhi("hEzZPhi", "Local longitudinal field for PHGarfield averaged over r;|z| [cm];#phi [rad]", m_nzObs, zcm.data(), m_nphiObs, phi_edges.data());
  hEr.GetXaxis()->SetTitle("r [cm]"); hEr.GetYaxis()->SetTitle("|z| [cm]"); hEr.GetZaxis()->SetTitle("E_{r} [V/m]");
  hEphi.GetXaxis()->SetTitle("r [cm]"); hEphi.GetYaxis()->SetTitle("|z| [cm]"); hEphi.GetZaxis()->SetTitle("E_{#phi} [V/m]");
  hEz.GetXaxis()->SetTitle("r [cm]"); hEz.GetYaxis()->SetTitle("|z| [cm]"); hEz.GetZaxis()->SetTitle("E_{z}^{local} [V/m]");
  hErRPhi.GetXaxis()->SetTitle("r [cm]"); hErRPhi.GetYaxis()->SetTitle("#phi [rad]"); hErRPhi.GetZaxis()->SetTitle("E_{r} [V/m]");
  hEphiRPhi.GetXaxis()->SetTitle("r [cm]"); hEphiRPhi.GetYaxis()->SetTitle("#phi [rad]"); hEphiRPhi.GetZaxis()->SetTitle("E_{#phi} [V/m]");
  hEzRPhi.GetXaxis()->SetTitle("r [cm]"); hEzRPhi.GetYaxis()->SetTitle("#phi [rad]"); hEzRPhi.GetZaxis()->SetTitle("E_{z}^{local} [V/m]");
  hErZPhi.GetXaxis()->SetTitle("|z| [cm]"); hErZPhi.GetYaxis()->SetTitle("#phi [rad]"); hErZPhi.GetZaxis()->SetTitle("E_{r} [V/m]");
  hEphiZPhi.GetXaxis()->SetTitle("|z| [cm]"); hEphiZPhi.GetYaxis()->SetTitle("#phi [rad]"); hEphiZPhi.GetZaxis()->SetTitle("E_{#phi} [V/m]");
  hEzZPhi.GetXaxis()->SetTitle("|z| [cm]"); hEzZPhi.GetYaxis()->SetTitle("#phi [rad]"); hEzZPhi.GetZaxis()->SetTitle("E_{z}^{local} [V/m]");

  for (unsigned int ir = r_begin; ir < r_end; ++ir)
  {
    for (unsigned int iz = 0; iz < m_nzObs; ++iz)
    {
      double ersum = 0.0;
      double ephisum = 0.0;
      double ezsum = 0.0;
      for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        ersum += er[idx];
        ephisum += ephi[idx];
        ezsum += ez[idx];
      }
      hEr.SetBinContent(ir + 1, iz + 1, ersum / static_cast<double>(m_nphiObs));
      hEphi.SetBinContent(ir + 1, iz + 1, ephisum / static_cast<double>(m_nphiObs));
      hEz.SetBinContent(ir + 1, iz + 1, ezsum / static_cast<double>(m_nphiObs));
    }

    for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
    {
      double ersum = 0.0;
      double ephisum = 0.0;
      double ezsum = 0.0;
      for (unsigned int iz = 0; iz < m_nzObs; ++iz)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        ersum += er[idx];
        ephisum += ephi[idx];
        ezsum += ez[idx];
      }
      hErRPhi.SetBinContent(ir + 1, ip + 1, ersum / static_cast<double>(m_nzObs));
      hEphiRPhi.SetBinContent(ir + 1, ip + 1, ephisum / static_cast<double>(m_nzObs));
      hEzRPhi.SetBinContent(ir + 1, ip + 1, ezsum / static_cast<double>(m_nzObs));
    }
  }

  const unsigned int nr_filled = r_end - r_begin;
  if (nr_filled > 0)
  {
    for (unsigned int iz = 0; iz < m_nzObs; ++iz)
    {
      for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      {
        double ersum = 0.0;
        double ephisum = 0.0;
        double ezsum = 0.0;
        for (unsigned int ir = r_begin; ir < r_end; ++ir)
        {
          const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
          ersum += er[idx];
          ephisum += ephi[idx];
          ezsum += ez[idx];
        }
        hErZPhi.SetBinContent(iz + 1, ip + 1, ersum / static_cast<double>(nr_filled));
        hEphiZPhi.SetBinContent(iz + 1, ip + 1, ephisum / static_cast<double>(nr_filled));
        hEzZPhi.SetBinContent(iz + 1, ip + 1, ezsum / static_cast<double>(nr_filled));
      }
    }
  }

  hEr.Write();
  hEphi.Write();
  hEz.Write();
  hErRPhi.Write();
  hEphiRPhi.Write();
  hEzRPhi.Write();
  hErZPhi.Write();
  hEphiZPhi.Write();
  hEzZPhi.Write();
  file.cd();
  const std::string meta = std::format("{{\"format\":\"PHGarfield axisymmetric correction map\",\"histograms\":{{\"Er\":\"QA/hErDefault\",\"Ephi\":\"QA/hEphiDefault\",\"Ez\":\"QA/hEzDefault\",\"ErRPhi\":\"QA/hErRPhi\",\"EphiRPhi\":\"QA/hEphiRPhi\",\"EzRPhi\":\"QA/hEzRPhi\",\"ErZPhi\":\"QA/hErZPhi\",\"EphiZPhi\":\"QA/hEphiZPhi\",\"EzZPhi\":\"QA/hEzZPhi\"}},\"field_units\":\"V/m\",\"job_index\":{},\"n_jobs\":{},\"radial_bin_begin\":{},\"radial_bin_end\":{}}}", m_jobIndex, m_nJobs, r_begin, r_end);
  TNamed("metadata_json", meta.c_str()).Write();
  file.Close();
  std::cout << Name() << " wrote " << m_garfieldOutputFile << std::endl;
}

void PHGarfieldRossegger::writeField3DRootFile(const std::vector<double>& r_source_edges_m, const std::vector<double>& r_obs_edges_m, const std::vector<double>& phi_source_edges, const std::vector<double>& phi_obs_edges, const std::vector<double>& z_source_edges_m, const std::vector<double>& z_obs_edges_m, const std::vector<double>& rho, const std::vector<double>& potential, const std::vector<double>& er, const std::vector<double>& ephi, const std::vector<double>& ez, unsigned int r_begin, unsigned int r_end) const
{
  auto rscm = r_source_edges_m;
  auto rocm = r_obs_edges_m;
  auto zscm = z_source_edges_m;
  auto zocm = z_obs_edges_m;
  scale_to_cm(rscm); scale_to_cm(rocm); scale_to_cm(zscm); scale_to_cm(zocm);
  TFile file(m_field3DOutputFile.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + m_field3DOutputFile); }
  const unsigned int nr_source = r_source_edges_m.size() - 1;
  TH3D hRho("hRho", "Positive-ion density", nr_source, rscm.data(), m_nphiSource, phi_source_edges.data(), m_nzSource, zscm.data());
  TH3D hPhi("hPhi", "Space-charge potential", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEr("hEr", "Radial electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEphi("hEphi", "Azimuthal electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEz("hEz", "Longitudinal electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  for (unsigned int ir = 0; ir < nr_source; ++ir)
  {
    for (unsigned int ip = 0; ip < m_nphiSource; ++ip)
    {
      for (unsigned int iz = 0; iz < m_nzSource; ++iz)
      {
        const unsigned int idx = src_index(ir, ip, iz, m_nphiSource, m_nzSource);
        hRho.SetBinContent(ir + 1, ip + 1, iz + 1, rho[idx]);
      }
    }
  }
  for (unsigned int ir = r_begin; ir < r_end; ++ir)
    for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      for (unsigned int iz = 0; iz < m_nzObs; ++iz)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        hPhi.SetBinContent(ir + 1, ip + 1, iz + 1, potential[idx]);
        hEr.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx]);
        hEphi.SetBinContent(ir + 1, ip + 1, iz + 1, ephi[idx]);
        hEz.SetBinContent(ir + 1, ip + 1, iz + 1, ez[idx]);
      }
  hRho.Write(); hPhi.Write(); hEr.Write(); hEphi.Write(); hEz.Write();
  file.cd();
  TNamed("metadata_json", "Rossegger modal solution for sPHENIX TPC").Write();
  file.Close();
  std::cout << Name() << " wrote " << m_field3DOutputFile << std::endl;
}


void PHGarfieldRossegger::writeFrameBoundaryField3DRootFile(const std::vector<double>& r_source_edges_m,
                                                            const std::vector<double>& r_obs_edges_m,
                                                            const std::vector<double>& phi_source_edges,
                                                            const std::vector<double>& phi_obs_edges,
                                                            const std::vector<double>& z_obs_edges_m,
                                                            const FrameBoundaryPattern& frame_pattern,
                                                            const std::vector<double>& potential,
                                                            const std::vector<double>& er,
                                                            const std::vector<double>& ephi,
                                                            const std::vector<double>& ez,
                                                            unsigned int r_begin,
                                                            unsigned int r_end) const
{
  auto rscm = r_source_edges_m;
  auto rocm = r_obs_edges_m;
  auto zocm = z_obs_edges_m;
  scale_to_cm(rscm); scale_to_cm(rocm); scale_to_cm(zocm);
  TFile file(m_field3DOutputFile.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + m_field3DOutputFile); }

  const unsigned int nr_source = r_source_edges_m.size() - 1;
  TH2D hFrameGeometry("hFrameGeometry", "Frame geometry fraction at GEM;r [cm];#phi [rad]", nr_source, rscm.data(), m_nphiSource, phi_source_edges.data());
  TH2D hFrameWeight("hFrameWeight", "Frame boundary weighting at GEM;r [cm];#phi [rad]", nr_source, rscm.data(), m_nphiSource, phi_source_edges.data());
  TH2D hFrameBoundaryPotential("hFrameBoundaryPotential", "Effective frame boundary potential at GEM;r [cm];#phi [rad]", nr_source, rscm.data(), m_nphiSource, phi_source_edges.data());
  TH3D hPhi("hPhi", "Frame boundary potential solution;r [cm];#phi [rad];|z| [cm]", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEr("hEr", "Radial electric field from frame boundary;r [cm];#phi [rad];|z| [cm]", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEphi("hEphi", "Azimuthal electric field from frame boundary;r [cm];#phi [rad];|z| [cm]", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEz("hEz", "Longitudinal electric field from frame boundary;r [cm];#phi [rad];|z| [cm]", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());

  for (unsigned int ir = 0; ir < nr_source; ++ir)
  {
    for (unsigned int ip = 0; ip < m_nphiSource; ++ip)
    {
      const unsigned int idx = ir * m_nphiSource + ip;
      hFrameGeometry.SetBinContent(ir + 1, ip + 1, frame_pattern.geometry_fraction[idx]);
      hFrameWeight.SetBinContent(ir + 1, ip + 1, frame_pattern.weight[idx]);
      hFrameBoundaryPotential.SetBinContent(ir + 1, ip + 1, frame_pattern.boundary_potential[idx]);
    }
  }

  for (unsigned int ir = r_begin; ir < r_end; ++ir)
    for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      for (unsigned int iz = 0; iz < m_nzObs; ++iz)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        hPhi.SetBinContent(ir + 1, ip + 1, iz + 1, potential[idx]);
        hEr.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx]);
        hEphi.SetBinContent(ir + 1, ip + 1, iz + 1, ephi[idx]);
        hEz.SetBinContent(ir + 1, ip + 1, iz + 1, ez[idx]);
      }

  hFrameGeometry.Write();
  hFrameWeight.Write();
  hFrameBoundaryPotential.Write();
  hPhi.Write(); hEr.Write(); hEphi.Write(); hEz.Write();
  TNamed("metadata_json", std::format("{{\"format\":\"Rossegger frame boundary solution\",\"frame_boundary_potential_V\":{},\"frame_charge_weighting\":\"{}\"}}", m_frameBoundaryPotential, m_frameChargeWeighting == FrameChargeWeighting::EqualChargePerPiece ? "EqualChargePerPiece" : "ProportionalToArea").c_str()).Write();
  file.Close();
  std::cout << Name() << " wrote " << m_field3DOutputFile << std::endl;
}

void PHGarfieldRossegger::writePHGarfieldField3DRootFile(const std::vector<double>& r_obs_edges_m,
                                                         const std::vector<double>& phi_obs_edges,
                                                         const std::vector<double>& z_obs_edges_m,
                                                         const std::vector<double>& potential,
							 const std::vector<double>& er,
                                                         const std::vector<double>& ephi,
                                                         const std::vector<double>& ez,
                                                         unsigned int side,
                                                         unsigned int r_begin,
                                                         unsigned int r_end) const
{
  auto rcm = r_obs_edges_m;
  auto zcm = z_obs_edges_m;
  scale_to_cm(rcm);
  scale_to_cm(zcm);

  const std::string& filename = m_phGarfieldField3DOutputFiles[side];
  TFile file(filename.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + filename); }
  auto* dir = file.mkdir("Field3D");
  if (!dir) { throw std::runtime_error("Could not create Field3D in " + filename); }
  dir->cd();

  TH3D hPhi("hPhi",std::format("TPC electric potential, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(),m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  TH3D hEr("hEr",std::format("TPC radial electric-field correction, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(),m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  TH3D hEphi("hEphi",std::format("TPC azimuthal electric-field correction, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(),m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  TH3D hEx("hEx", std::format("Local TPC E_{{x}} correction, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(), m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  TH3D hEy("hEy", std::format("Local TPC E_{{y}} correction, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(), m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  TH3D hEz("hEz", std::format("TPC E_{{z}} correction in +|z| solver coordinate, side {};r [cm];#phi [rad];|z| [cm]", side).c_str(), m_nrObs, rcm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zcm.data());
  hEx.GetXaxis()->SetTitle("r [cm]"); hEx.GetYaxis()->SetTitle("#phi [rad]"); hEx.GetZaxis()->SetTitle("|z| [cm]");
  hEy.GetXaxis()->SetTitle("r [cm]"); hEy.GetYaxis()->SetTitle("#phi [rad]"); hEy.GetZaxis()->SetTitle("|z| [cm]");
  hEz.GetXaxis()->SetTitle("r [cm]"); hEz.GetYaxis()->SetTitle("#phi [rad]"); hEz.GetZaxis()->SetTitle("|z| [cm]");

  const auto phi_cent = centers(phi_obs_edges);
  for (unsigned int ir = r_begin; ir < r_end; ++ir)
  {
    for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
    {
      const double cphi = std::cos(phi_cent[ip]);
      const double sphi = std::sin(phi_cent[ip]);
      for (unsigned int iz = 0; iz < m_nzObs; ++iz)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
	hPhi.SetBinContent(ir + 1, ip + 1, iz + 1, potential[idx]);
	hEr.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx]);
	hEphi.SetBinContent(ir + 1, ip + 1, iz + 1, ephi[idx]);
        hEx.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx] * cphi - ephi[idx] * sphi);
        hEy.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx] * sphi + ephi[idx] * cphi);
        hEz.SetBinContent(ir + 1, ip + 1, iz + 1, ez[idx]);
      }
    }
  }

  hPhi.SetEntries(static_cast<double>(m_nrObs) * static_cast<double>(m_nphiObs) * static_cast<double>(m_nzObs));
  hEr.SetEntries(hPhi.GetEntries());
  hEphi.SetEntries(hPhi.GetEntries());
  hEx.SetEntries(static_cast<double>(m_nrObs) * static_cast<double>(m_nphiObs) * static_cast<double>(m_nzObs));
  hEy.SetEntries(hEx.GetEntries());
  hEz.SetEntries(hEx.GetEntries());
  hPhi.Write("hPhi", TObject::kOverwrite);
  hEr.Write("hEr", TObject::kOverwrite);
  hEphi.Write("hEphi", TObject::kOverwrite);
  hEx.Write("hEx", TObject::kOverwrite);
  hEy.Write("hEy", TObject::kOverwrite);
  hEz.Write("hEz", TObject::kOverwrite);

  file.cd();
  TNamed("field_map_format", "PHGarfield Field3D v1").Write();
  TNamed("coordinate_convention", std::format("local TPC Cartesian Ex,Ey,Ez for side {}; axes=(r_cm,phi_rad,abs_z_cm); phi in [-pi,pi], periodic; histogram contents V/m; Ez is derivative along the +|z| solver coordinate", side).c_str()).Write();
  const std::string meta = std::format("{{\"format\":\"PHGarfield Field3D v1\",\"coordinate_system\":\"local TPC\",\"side\":{},\"gain_file\":\"{}\",\"gain_histogram\":\"{}\",\"gain_application\":\"rho_corrected = rho_base / gain\",\"field_components\":[\"Ex\",\"Ey\",\"Ez\"],\"field_units\":\"V/m\",\"r_axis_units\":\"cm\",\"phi_axis_units\":\"rad\",\"z_axis_units\":\"cm\",\"z_axis_definition\":\"abs(z); Ez stored in +|z| solver coordinate\",\"phi_periodic\":true,\"phi_range\":[{},{}],\"shape\":[{},{},{}],\"source_phi_cells\":{},\"source_z_cells\":{},\"N_RADIAL_MODES\":{},\"N_LONGITUDINAL_MODES\":{},\"M_PHI_MAX\":{},\"K_EFF\":{}}}", side, m_gainMapFile, m_gainHistograms[side], phi_obs_edges.front(), phi_obs_edges.back(), m_nrObs, m_nphiObs, m_nzObs, m_nphiSource, m_nzSource, m_nRadialModes, m_nLongitudinalModes, m_mPhiMax, m_kEff);
  TNamed("metadata_json", meta.c_str()).Write();
  file.Close();
  std::cout << Name() << " wrote " << filename << std::endl;
}

bool PHGarfieldRossegger::verifyOutput() const
{
  TFile file(m_garfieldOutputFile.c_str(), "READ");
  if (!file.IsOpen() || file.IsZombie()) { return false; }
  return file.Get("QA/hErDefault") && file.Get("QA/hEphiDefault") && file.Get("QA/hEzDefault") &&
         file.Get("QA/hErRPhi") && file.Get("QA/hEphiRPhi") && file.Get("QA/hEzRPhi") &&
         file.Get("QA/hErZPhi") && file.Get("QA/hEphiZPhi") && file.Get("QA/hEzZPhi") &&
         file.Get("metadata_json");
}

#include "PHGarfieldForwardMap.h"

#include <TAxis.h>
#include <TFile.h>
#include <TH3.h>
#include <TParameter.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

namespace
{
  const std::string MetadataPath = "GarfieldFullTransport/Metadata";
  void require(bool valid, const std::string& message)
  {
    if (!valid)
    {
      throw std::runtime_error(message);
    }
  }
  template <class T>
  T parameter(TDirectory& directory, const char* name)
  {
    const auto* value = dynamic_cast<TParameter<T>*>(directory.Get(name));
    require(value != nullptr, std::string("missing/incorrect metadata ") + name);
    return value->GetVal();
  }
  bool close_value(const double first, const double second)
  {
    return std::isfinite(first) && std::isfinite(second) &&
           std::abs(first - second) <= 1.e-12 * std::max({1., std::abs(first), std::abs(second)});
  }

  void check_axis(const TAxis& axis, int bins, double minimum, double maximum)
  {
    require(axis.GetNbins() == bins, "histogram dimensions disagree with physical metadata");
    for (int edge = 1; edge <= bins + 1; ++edge)
    {
      const double expected = minimum + (maximum - minimum) * (edge - 1) / bins;
      require(close_value(axis.GetBinLowEdge(edge), expected), "histogram axis disagrees with physical metadata");
    }
  }

  void validate_merged_map(TFile& file, std::ostream& error)
  {
    auto* metadata = file.GetDirectory(MetadataPath.c_str());
    require(metadata != nullptr, "merged metadata missing");
    PHGarfieldForwardMap::Config geometry;
    geometry.output_filename = file.GetName();
    const auto positive = [&](const char* name)
    {
      const int value = parameter<int>(*metadata, name);
      require(value > 0, std::string("nonpositive geometry ") + name);
      return static_cast<unsigned int>(value);
    };
    geometry.phi_nodes = positive("phi_nodes");
    geometry.radial_nodes = positive("radial_nodes");
    geometry.longitudinal_nodes = positive("longitudinal_nodes");
    geometry.inner_radius_cm = parameter<double>(*metadata, "gas_inner_radius_cm");
    geometry.outer_radius_cm = parameter<double>(*metadata, "gas_outer_radius_cm");
    geometry.central_membrane_half_width_cm = parameter<double>(*metadata, "central_membrane_half_width_cm");
    geometry.readout_abs_z_cm = parameter<double>(*metadata, "readout_abs_z_cm");
    geometry.integration_step_ns = parameter<double>(*metadata, "integration_step_ns");
    geometry.nominal_phg4_drift_velocity_cm_per_ns = parameter<double>(*metadata, "nominal_phg4_drift_velocity_cm_per_ns");
    geometry.maximum_nominal_time_mismatch = parameter<double>(*metadata, "maximum_nominal_time_mismatch");
    geometry.probe_cell_centers = parameter<int>(*metadata, "cell_center_probes") == 1;
    require(geometry.probe_cell_centers, "runtime maps require interior reach probes");
    const int version = parameter<int>(*metadata, "schema_version");
    if (version >= PHGarfieldForwardMap::NumericalSchemaVersion)
    {
      geometry.interpolation_tolerance_um = parameter<double>(*metadata, "interpolation_tolerance_um");
      require(geometry.interpolation_tolerance_um > 0., "schema3 requires a positive interpolation tolerance");
    }
    require(PHGarfieldForwardMap::ValidateConfig(geometry, error), "invalid merged geometry/numerical configuration");
    const double maximum_time = parameter<double>(*metadata, "maximum_drift_time_ns");
    require(std::isfinite(maximum_time) && maximum_time >= 0., "invalid merged maximum drift time");
    double observed_maximum_time = 0.;
    const double pi = std::acos(-1.);
    const auto get_histogram = [&](const std::string& stem, const std::string& side) -> const TH3&
    {
      const auto name = "GarfieldFullTransport/" + stem + "_" + side;
      const auto* histogram = dynamic_cast<const TH3*>(file.Get(name.c_str()));
      require(histogram != nullptr, "missing or incorrect histogram type: " + name);
      return *histogram;
    };
    for (const std::string side : {"negz", "posz"})
    {
      const double minimum_z = side == "negz" ? -geometry.readout_abs_z_cm : geometry.central_membrane_half_width_cm;
      const double maximum_z = side == "negz" ? -geometry.central_membrane_half_width_cm : geometry.readout_abs_z_cm;
      const auto& radial = get_histogram("hDeltaR", side);
      const auto& tangential = get_histogram("hDeltaRPhi", side);
      const auto& duration = get_histogram("hDriftTime", side);
      const auto& reaches = get_histogram("hReachesReadoutNode", side);
      const auto& cells = get_histogram("hReachCellClass", side);
      const double phi_step = 2. * pi / geometry.phi_nodes;
      const double radius_step = (geometry.outer_radius_cm - geometry.inner_radius_cm) / (geometry.radial_nodes - 1);
      const double z_step = (maximum_z - minimum_z) / (geometry.longitudinal_nodes - 1);
      for (const auto* histogram : {&radial, &tangential, &duration, &reaches})
      {
        check_axis(*histogram->GetXaxis(), geometry.phi_nodes + 2, -1.5 * phi_step, 2. * pi + .5 * phi_step);
        check_axis(*histogram->GetYaxis(), geometry.radial_nodes + 2, geometry.inner_radius_cm - 1.5 * radius_step, geometry.outer_radius_cm + 1.5 * radius_step);
        check_axis(*histogram->GetZaxis(), geometry.longitudinal_nodes + 2, minimum_z - 1.5 * z_step, maximum_z + 1.5 * z_step);
        // Periodic phi guards and nearest physical radial/z support must be
        // copied exactly, including combined corners of the guard region.
        for (int x = 1; x <= histogram->GetNbinsX(); ++x)
        {
          for (int y = 1; y <= histogram->GetNbinsY(); ++y)
          {
            for (int z = 1; z <= histogram->GetNbinsZ(); ++z)
            {
              const int source_x = x == 1 ? histogram->GetNbinsX() - 1 : (x == histogram->GetNbinsX() ? 2 : x);
              const int source_y = std::clamp(y, 2, histogram->GetNbinsY() - 1);
              const int source_z = std::clamp(z, 2, histogram->GetNbinsZ() - 1);
              require(histogram->GetBinContent(x, y, z) == histogram->GetBinContent(source_x, source_y, source_z), "invalid periodic or copied guard node");
            }
          }
        }
      }
      const auto check_cell_axes = [&](const TH3& histogram)
      {
        check_axis(*histogram.GetXaxis(), geometry.phi_nodes, 0., 2. * pi);
        check_axis(*histogram.GetYaxis(), geometry.radial_nodes - 1, geometry.inner_radius_cm, geometry.outer_radius_cm);
        check_axis(*histogram.GetZaxis(), geometry.longitudinal_nodes - 1, minimum_z, maximum_z);
      };
      check_cell_axes(cells);
      const TH3* transverse_error = nullptr;
      const TH3* time_error = nullptr;
      if (version >= PHGarfieldForwardMap::NumericalSchemaVersion)
      {
        transverse_error = &get_histogram("hInterpolationTransverseError", side);
        time_error = &get_histogram("hInterpolationTimeEquivalentError", side);
        check_cell_axes(*transverse_error);
        check_cell_axes(*time_error);
      }
      for (int x = 2; x < radial.GetNbinsX(); ++x)
      {
        for (int y = 2; y < radial.GetNbinsY(); ++y)
        {
          for (int z = 2; z < radial.GetNbinsZ(); ++z)
          {
            const double reach = reaches.GetBinContent(x, y, z);
            const double dr = radial.GetBinContent(x, y, z);
            const double rdphi = tangential.GetBinContent(x, y, z);
            const double time = duration.GetBinContent(x, y, z);
            require((reach == 0. || reach == 1.) && std::isfinite(dr) && std::isfinite(rdphi) && std::isfinite(time) && time >= 0., "invalid node reach or continuous transport value");
            require(reach != 0. || (dr == 0. && rdphi == 0. && time == 0.), "lost node has nonzero transport");
            if (reach == 1.)
            {
              const double radius = radial.GetYaxis()->GetBinCenter(y);
              require(radius + dr >= geometry.inner_radius_cm - 1.e-9 && radius + dr <= geometry.outer_radius_cm + 1.e-9 && std::abs(rdphi / radius) <= pi + 1.e-9 && time <= maximum_time + 1.e-9, "reaching node violates physical endpoint/time semantics");
              observed_maximum_time = std::max(observed_maximum_time, time);
            }
            const int readout_z = side == "negz" ? 2 : radial.GetNbinsZ() - 1;
            require(z != readout_z || (reach == 1. && std::abs(dr) <= 1.e-12 && std::abs(rdphi) <= 1.e-12 && time <= 1.e-12), "readout node is not zero transport");
          }
        }
      }
      for (int x = 1; x <= cells.GetNbinsX(); ++x)
      {
        for (int y = 1; y <= cells.GetNbinsY(); ++y)
        {
          for (int z = 1; z <= cells.GetNbinsZ(); ++z)
          {
            const double code = cells.GetBinContent(x, y, z);
            require(code == 0. || code == 1. || code == 2. || (code == 3. && version >= PHGarfieldForwardMap::NumericalSchemaVersion), "invalid reach-cell code");
            unsigned int reaching_corners = 0;
            double minimum_phi = std::numeric_limits<double>::infinity();
            double maximum_phi = -minimum_phi;
            for (int cx : {x + 1, x + 2})
            {
              for (int cy : {y + 1, y + 2})
              {
                for (int cz : {z + 1, z + 2})
                {
                  if (reaches.GetBinContent(cx, cy, cz) == 1.)
                  {
                    ++reaching_corners;
                    const double phi = tangential.GetBinContent(cx, cy, cz) / reaches.GetYaxis()->GetBinCenter(cy);
                    minimum_phi = std::min(minimum_phi, phi);
                    maximum_phi = std::max(maximum_phi, phi);
                  }
                }
              }
            }
            require((code != 0. || reaching_corners == 0) && ((code != 2. && code != 3.) || (reaching_corners == 8 && maximum_phi - minimum_phi <= pi)), "reach-cell classification disagrees with its support corners");
            if (transverse_error)
            {
              const double transverse = transverse_error->GetBinContent(x, y, z);
              const double time = time_error->GetBinContent(x, y, z);
              require(std::isfinite(transverse) && transverse >= 0. && std::isfinite(time) && time >= 0., "invalid interpolation-error diagnostic");
              require((code != 2. && code != 3.) || ((std::max(transverse, time) > geometry.interpolation_tolerance_um) == (code == 3.)), "numerical-mask class disagrees with sampled error");
            }
          }
        }
      }
    }
    require(close_value(observed_maximum_time, maximum_time), "maximum drift-time metadata disagrees with source nodes");
  }

}  // namespace

bool PHGarfieldForwardMap::ValidateRegular(TFile& file, std::ostream& error)
{
  try
  {
    validate_merged_map(file, error);
    return true;
  }
  catch (const std::exception& exception)
  {
    error << "Invalid regular forward LUT: " << exception.what() << '\n';
    return false;
  }
}

#include "PHGarfieldForwardMap.h"
#include "PHGarfieldFingerprint.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TH3I.h>
#include <TMD5.h>
#include <TNamed.h>
#include <TParameter.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <regex>
#include <sstream>
#include <string>
#include <system_error>

namespace
{
  constexpr double TwoPi = 2.0 * std::numbers::pi;

  struct UniformNodeAxis
  {
    unsigned int physical_nodes{0};
    double physical_min{0.0};
    double physical_max{0.0};
    double spacing{0.0};
    double histogram_min{0.0};
    double histogram_max{0.0};
  };

  UniformNodeAxis make_node_axis(
      const unsigned int nodes,
      const double minimum,
      const double maximum)
  {
    UniformNodeAxis axis;
    axis.physical_nodes = nodes;
    axis.physical_min = minimum;
    axis.physical_max = maximum;
    axis.spacing = (maximum - minimum) / static_cast<double>(nodes - 1);
    axis.histogram_min = minimum - 1.5 * axis.spacing;
    axis.histogram_max = maximum + 1.5 * axis.spacing;
    return axis;
  }

  UniformNodeAxis make_phi_axis(const unsigned int nodes)
  {
    const double spacing = TwoPi / static_cast<double>(nodes);
    UniformNodeAxis axis;
    axis.physical_nodes = nodes;
    axis.physical_min = 0.0;
    axis.physical_max = TwoPi - spacing;
    axis.spacing = spacing;
    axis.histogram_min = -1.5 * spacing;
    axis.histogram_max = TwoPi + 0.5 * spacing;
    return axis;
  }

  double physical_node_value(
      const UniformNodeAxis& axis,
      const unsigned int index)
  {
    if (index == 0)
    {
      return axis.physical_min;
    }
    if (index + 1 == axis.physical_nodes)
    {
      return axis.physical_max;
    }
    return axis.physical_min + static_cast<double>(index) * axis.spacing;
  }

  template <class Histogram>
  void fill_guard_nodes(Histogram& histogram,
                        const unsigned int phi_nodes,
                        const unsigned int radial_nodes,
                        const unsigned int longitudinal_nodes)
  {
    const int nx = static_cast<int>(phi_nodes) + 2;
    const int ny = static_cast<int>(radial_nodes) + 2;
    const int nz = static_cast<int>(longitudinal_nodes) + 2;

    for (int ix = 1; ix <= nx; ++ix)
    {
      const int source_x =
          ix == 1
              ? static_cast<int>(phi_nodes) + 1
              : (ix == nx ? 2 : ix);
      for (int iy = 1; iy <= ny; ++iy)
      {
        const int source_y =
            iy == 1
                ? 2
                : (iy == ny ? static_cast<int>(radial_nodes) + 1 : iy);
        for (int iz = 1; iz <= nz; ++iz)
        {
          if (ix != 1 && ix != nx &&
              iy != 1 && iy != ny &&
              iz != 1 && iz != nz)
          {
            continue;
          }

          const int source_z =
              iz == 1
                  ? 2
                  : (iz == nz
                         ? static_cast<int>(longitudinal_nodes) + 1
                         : iz);
          histogram.SetBinContent(
              ix,
              iy,
              iz,
              histogram.GetBinContent(source_x, source_y, source_z));
        }
      }
    }
  }

  bool physical_loss(const PHGarfield::ForwardDriftResult& result)
  {
    return result.status == PHGarfield::ForwardDriftStatus::LeftDriftVolume;
  }

  bool run_drift(const PHGarfield& garfield,
                 const double r_cm,
                 const double phi_rad,
                 const double z_cm,
                 const PHGarfieldForwardMap::Config& config,
                 PHGarfieldForwardMap::NodeValue& value,
                 std::ostream& error)
  {
    PHGarfield::ForwardDriftResult result;
    const bool ok = garfield.ForwardDrift(
        r_cm * std::cos(phi_rad),
        r_cm * std::sin(phi_rad),
        z_cm,
        config.readout_abs_z_cm,
        config.integration_step_ns,
        result,
        true);

    if (ok && result.reached_readout)
    {
      if (!PHGarfieldForwardMap::EncodeNode(r_cm, phi_rad, result, value))
      {
        error << "Failed to encode a reaching PHGarfield trajectory at (r,phi,z)=("
              << r_cm << ", " << phi_rad << ", " << z_cm << ")\n";
        return false;
      }
      return true;
    }

    if (physical_loss(result))
    {
      value = PHGarfieldForwardMap::NodeValue{};
      return true;
    }

    error << "Fatal PHGarfield trajectory failure at (r,phi,z)=("
          << r_cm << ", " << phi_rad << ", " << z_cm << "): "
          << PHGarfield::ForwardDriftStatusName(result.status) << "\n";
    return false;
  }

  std::string regular_file_md5(const std::filesystem::path& filename)
  {
    std::unique_ptr<TMD5> checksum(
        TMD5::FileChecksum(filename.string().c_str()));
    return checksum ? checksum->AsString() : std::string{};
  }

  template <class T>
  void write_parameter(const char* name, const T value)
  {
    TParameter<T> parameter(name, value);
    parameter.Write();
  }

  void write_named(const char* name, const std::string& value)
  {
    TNamed named(name, value.c_str());
    named.Write();
  }

}  // namespace

std::string PHGarfieldForwardMap::PayloadMd5(const std::string& filename)
{
  try
  {
    if (filename.empty())
    {
      return {};
    }
    const std::filesystem::path path(filename);
    std::error_code error;
    if (std::filesystem::is_regular_file(path, error) && !error)
    {
      return regular_file_md5(path);
    }
    if (!std::filesystem::is_directory(path, error) || error)
    {
      return {};
    }

    // Match PHGarfield::InitializeGas: only ordered MERGED_E###.gas payloads
    // contribute to a directory-based gas table.
    const std::regex gas_pattern(R"(^MERGED_E([0-9]{3})\.gas$)");
    std::map<unsigned int, std::filesystem::path> gas_files;
    for (const auto& entry : std::filesystem::directory_iterator(path))
    {
      std::smatch match;
      const std::string basename = entry.path().filename().string();
      if (entry.is_regular_file() &&
          std::regex_match(basename, match, gas_pattern))
      {
        gas_files[std::stoul(match[1].str())] = entry.path();
      }
    }
    if (gas_files.empty())
    {
      return {};
    }

    std::ostringstream manifest;
    for (const auto& [field_index, gas_path] : gas_files)
    {
      const std::string digest = regular_file_md5(gas_path);
      if (digest.empty())
      {
        return {};
      }
      manifest << field_index << ":" << gas_path.filename().string()
               << ":" << digest << "\n";
    }
    const std::string text = manifest.str();
    TMD5 checksum;
    checksum.Update(
        reinterpret_cast<const unsigned char*>(text.data()),
        static_cast<unsigned int>(text.size()));
    checksum.Final();
    return checksum.AsString();
  }
  catch (const std::filesystem::filesystem_error&)
  {
    return {};
  }
}

std::string PHGarfieldForwardMap::Fingerprint(const std::string& text)
{
  PHGarfieldFingerprint::Fnv1a64 hash;
  hash.Bytes(text);
  return hash.Hex();
}

bool PHGarfieldForwardMap::EncodeNode(
    const double r0_cm,
    const double phi0_rad,
    const PHGarfield::ForwardDriftResult& result,
    NodeValue& value)
{
  value = NodeValue{};
  if (!result.reached_readout ||
      result.status != PHGarfield::ForwardDriftStatus::ReachedReadout ||
      !std::isfinite(r0_cm) ||
      !std::isfinite(phi0_rad) ||
      !std::isfinite(result.x_cm) ||
      !std::isfinite(result.y_cm) ||
      !std::isfinite(result.drift_time_ns) ||
      r0_cm <= 0.0 ||
      result.drift_time_ns < 0.0)
  {
    return false;
  }

  const double endpoint_radius_cm = std::hypot(result.x_cm, result.y_cm);
  const double endpoint_phi_rad = std::atan2(result.y_cm, result.x_cm);
  const double delta_phi_rad =
      std::remainder(endpoint_phi_rad - phi0_rad, TwoPi);

  value.reaches_readout = true;
  value.delta_r_cm = endpoint_radius_cm - r0_cm;
  value.delta_rphi_cm = r0_cm * delta_phi_rad;
  value.drift_time_ns = result.drift_time_ns;

  return std::isfinite(value.delta_r_cm) &&
         std::isfinite(value.delta_rphi_cm);
}

bool PHGarfieldForwardMap::ValidateConfig(
    const Config& config,
    std::ostream& error)
{
  bool valid = true;
  const auto require = [&valid, &error](const bool condition,
                                        const char* message)
  {
    if (!condition)
    {
      error << message << "\n";
      valid = false;
    }
  };

  require(!config.output_filename.empty(), "Output filename is empty");
  require(std::isfinite(config.interpolation_tolerance_um) && config.interpolation_tolerance_um >= 0., "Interpolation tolerance must be finite/nonnegative");
  require(config.interpolation_tolerance_um == 0. || config.probe_cell_centers, "Numerical mask requires interior probes");
  require(config.phi_nodes >= 3, "At least three phi nodes are required");
  require(config.radial_nodes >= 2, "At least two radial nodes are required");
  require(config.phi_chunk_count >= 1 && config.phi_chunk_count <= config.phi_nodes,
          "Phi chunk count must be within the number of physical phi nodes");
  require(config.phi_chunk_index < config.phi_chunk_count, "Phi chunk index is out of range");
  require(config.longitudinal_nodes >= 2,
          "At least two longitudinal nodes are required");
  constexpr unsigned int maximum_nodes_per_axis = 100000;
  require(
      config.phi_nodes <= maximum_nodes_per_axis &&
          config.radial_nodes <= maximum_nodes_per_axis &&
          config.longitudinal_nodes <= maximum_nodes_per_axis,
      "Grid node count exceeds the supported per-axis limit");
  constexpr std::uint64_t maximum_source_nodes = 10000000ULL;
  const std::uint64_t source_nodes =
      2ULL * static_cast<std::uint64_t>(config.phi_nodes) *
      static_cast<std::uint64_t>(config.radial_nodes) *
      static_cast<std::uint64_t>(config.longitudinal_nodes);
  require(
      source_nodes <= maximum_source_nodes,
      "Grid exceeds the supported total source-node limit");
  require(std::isfinite(config.inner_radius_cm) &&
              config.inner_radius_cm > 0.0,
          "Inner radius must be finite and positive");
  require(std::isfinite(config.outer_radius_cm) &&
              config.outer_radius_cm > config.inner_radius_cm,
          "Outer radius must exceed the inner radius");
  require(std::isfinite(config.central_membrane_half_width_cm) &&
              config.central_membrane_half_width_cm > 0.0,
          "Central-membrane half-width must be finite and positive");
  require(std::isfinite(config.readout_abs_z_cm) &&
              config.readout_abs_z_cm >
                  config.central_membrane_half_width_cm,
          "Readout |z| must exceed the central-membrane half-width");
  require(std::isfinite(config.integration_step_ns) &&
              config.integration_step_ns > 0.0,
          "Euler integration step must be finite and positive");
  require(
      std::isfinite(config.nominal_phg4_drift_velocity_cm_per_ns) &&
          config.nominal_phg4_drift_velocity_cm_per_ns > 0.0,
      "Nominal PHG4 drift velocity must be finite and positive");
  require(
      std::isfinite(config.maximum_nominal_time_mismatch) &&
          config.maximum_nominal_time_mismatch >= 0.0,
      "Maximum nominal-time mismatch must be finite and nonnegative");
  require(std::isfinite(config.tpc_place_x_cm) &&
              std::isfinite(config.tpc_place_y_cm) &&
              std::isfinite(config.tpc_place_z_cm) &&
              std::isfinite(config.tpc_rotation_x_rad) &&
              std::isfinite(config.tpc_rotation_y_rad) &&
              std::isfinite(config.tpc_rotation_z_rad),
          "TPC placement and rotation must be finite");
  return valid;
}

std::string PHGarfieldForwardMap::BuildCanonicalConfiguration(
    const Config& config,
    const Provenance& provenance,
    const int schema_version)
{
  std::ostringstream canonical;
  const auto append_string =
      [&canonical](const char* key, const std::string& value)
  {
    canonical << ";" << key << "[" << value.size() << "]=" << value;
  };
  canonical << std::setprecision(17)
            << "schema=" << SchemaName
            << ";version=" << schema_version
            << ";contract=" << TransportContract
            << ";phi_nodes=" << config.phi_nodes
            << ";radial_nodes=" << config.radial_nodes
            << ";longitudinal_nodes=" << config.longitudinal_nodes
            << ";inner_radius_cm=" << config.inner_radius_cm
            << ";outer_radius_cm=" << config.outer_radius_cm
            << ";central_membrane_half_width_cm="
            << config.central_membrane_half_width_cm
            << ";readout_abs_z_cm=" << config.readout_abs_z_cm
            << ";integration_step_ns=" << config.integration_step_ns
            << ";nominal_phg4_drift_velocity_cm_per_ns="
            << config.nominal_phg4_drift_velocity_cm_per_ns
            << ";maximum_nominal_time_mismatch="
            << config.maximum_nominal_time_mismatch
            << ";tpc_place_cm="
            << config.tpc_place_x_cm << ","
            << config.tpc_place_y_cm << ","
            << config.tpc_place_z_cm
            << ";tpc_rotation_rad="
            << config.tpc_rotation_x_rad << ","
            << config.tpc_rotation_y_rad << ","
            << config.tpc_rotation_z_rad
            << ";nominal_electric_field_vcm="
            << provenance.nominal_electric_field_vcm
            << ";space_charge_scale_side0="
            << provenance.space_charge_scale_side0
            << ";space_charge_scale_side1="
            << provenance.space_charge_scale_side1
            << ";zero_magnetic_field="
            << provenance.zero_magnetic_field
            << ";magnetic_field_scale="
            << provenance.magnetic_field_scale
            << ";magnetic_field_domain_policy="
            << MagneticFieldDomainPolicy
            << ";maximum_integration_steps="
            << provenance.maximum_integration_steps
            << ";cell_center_probes="
            << config.probe_cell_centers
            << ";nominal_time_reference_policy="
            << NominalTimeReferencePolicy
            << ";scale_zero_reference_time_side0_ns="
            << provenance.scale_zero_reference_time_side0_ns
            << ";scale_zero_reference_time_side1_ns="
            << provenance.scale_zero_reference_time_side1_ns
            << ";scale_zero_maximum_time_mismatch="
            << provenance.scale_zero_maximum_time_mismatch
            << ";cdb_timestamp=" << config.cdb_timestamp;
  append_string("electric_field_map", provenance.electric_field_map);
  append_string("magnetic_field_map", provenance.magnetic_field_map);
  append_string("gas_file", provenance.gas_file);
  append_string(
      "electric_field_md5",
      provenance.electric_field_map_md5);
  append_string(
      "magnetic_field_md5",
      provenance.magnetic_field_map_md5);
  append_string("gas_file_md5", provenance.gas_file_md5);
  append_string("reach_probe_scheme", provenance.reach_probe_scheme);
  append_string("cdb_global_tag", config.cdb_global_tag);
  append_string("configuration_label", config.configuration_label);
  if (schema_version >= ExtendedSchemaVersion)
  {
    append_string("field_model", provenance.field_model);
    if (provenance.field_model == "per_side_3d")
    {
      append_string("space_charge_z_coordinate_policy", SpaceChargeZCoordinatePolicy);
    }
    append_string("frame_field_map", provenance.frame_field_map);
    append_string("frame_field_map_md5", provenance.frame_field_map_md5);
    for (std::size_t side = 0; side < 2; ++side)
    {
      const std::string side_tag = "_side" + std::to_string(side);
      append_string(("electric_field_map_3d" + side_tag).c_str(), provenance.electric_field_map_3d[side]);
      append_string(("electric_field_map_3d_md5" + side_tag).c_str(), provenance.electric_field_map_3d_md5[side]);
      append_string(("frame_field_map_3d" + side_tag).c_str(), provenance.frame_field_map_3d[side]);
      append_string(("frame_field_map_3d_md5" + side_tag).c_str(), provenance.frame_field_map_3d_md5[side]);
      canonical << ";frame_charge_scale" << side_tag << "=" << provenance.frame_charge_scale[side]
                << ";cm_field_factor" << side_tag << "=" << provenance.cm_field_factor[side];
    }
    canonical << ";use_ifc_voltage_distortion=" << provenance.use_ifc_voltage_distortion
              << ";use_ofc_voltage_distortion=" << provenance.use_ofc_voltage_distortion;
    for (std::size_t index = 0; index < 4; ++index)
    {
      canonical << ";field_cage_offset_v_" << index << "=" << provenance.field_cage_offsets_v[index];
    }
  }
  if (schema_version >= NumericalSchemaVersion)
  {
    canonical << ";interpolation_tolerance_um=" << config.interpolation_tolerance_um;
    append_string("interpolation_error_probe_scheme", schema_version >= 4 ? provenance.reach_probe_scheme : "center_and_face_centers_transverse_time_equivalent_v1");
  }
  return canonical.str();
}

bool PHGarfieldForwardMap::Prepare(const PHGarfield& garfield,
                                   const Config& config, Provenance& provenance, Summary& summary,
                                   std::ostream& log, std::ostream& error)
{
  summary = Summary{};
  if (!ValidateConfig(config, error))
  {
    return false;
  }
  if (!garfield.IsReady())
  {
    error << "PHGarfield is not initialized\n";
    return false;
  }
  const bool has_3d = !garfield.GetElectricFieldMap3D(0).empty() || !garfield.GetElectricFieldMap3D(1).empty();
  if (config.interpolation_tolerance_um > 0. && !has_3d)
  {
    error << "Numerical-mask schema requires the explicit per-side3D model\n";
    return false;
  }
  if (garfield.GetFrameElectricFieldMap3D(0).empty() != garfield.GetFrameElectricFieldMap3D(1).empty())
  {
    error << "Both side-specific frame payloads are required\n";
    return false;
  }
  if ((garfield.GetSpaceChargeScaleSide0() != 0.0 &&
       !garfield.HasElectricFieldCorrections() && !garfield.HasLoadedElectricFieldMap3D(0)) ||
      (garfield.GetSpaceChargeScaleSide1() != 0.0 &&
       !garfield.HasElectricFieldCorrections() && !garfield.HasLoadedElectricFieldMap3D(1)))
  {
    error
        << "Nonzero electric-field correction scale requires a loaded "
           "electric-field correction map\n";
    return false;
  }

  constexpr double geometry_tolerance_cm = 1.0e-9;
  if (std::abs(garfield.GetDriftInnerRadius() - config.inner_radius_cm) >
          geometry_tolerance_cm ||
      std::abs(garfield.GetDriftOuterRadius() - config.outer_radius_cm) >
          geometry_tolerance_cm ||
      std::abs(garfield.GetCentralMembraneHalfWidth() -
               config.central_membrane_half_width_cm) >
          geometry_tolerance_cm)
  {
    error << "Generator grid geometry does not match the initialized "
             "PHGarfield drift geometry\n";
    return false;
  }
  if (!garfield.ForwardElectricFieldMapsCover(
          config.inner_radius_cm, config.outer_radius_cm,
          config.readout_abs_z_cm - config.central_membrane_half_width_cm, error))
  {
    return false;
  }
  if (!garfield.HasKnownMagneticFieldDomain() ||
      !garfield.MagneticFieldMapCoversTpcDriftVolume(
          config.outer_radius_cm,
          config.readout_abs_z_cm))
  {
    error
        << "The magnetic-field map domain does not conservatively cover the "
           "full placed TPC drift volume\n";
    return false;
  }

  const double reference_radius_cm =
      0.5 * (config.inner_radius_cm + config.outer_radius_cm);
  const double nominal_maximum_drift_time_ns =
      (config.readout_abs_z_cm -
       config.central_membrane_half_width_cm) /
      config.nominal_phg4_drift_velocity_cm_per_ns;
  std::array<double*, 2> reference_time = {
      {&summary.scale_zero_reference_time_side0_ns,
       &summary.scale_zero_reference_time_side1_ns}};
  for (unsigned int side = 0; side < 2; ++side)
  {
    PHGarfield::ForwardDriftResult reference_result;
    const double z_cm =
        side == 0
            ? -config.central_membrane_half_width_cm
            : config.central_membrane_half_width_cm;
    if (!garfield.ForwardDrift(
            reference_radius_cm,
            0.0,
            z_cm,
            config.readout_abs_z_cm,
            config.integration_step_ns,
            reference_result,
            false) ||
        !reference_result.reached_readout)
    {
      error
          << "Scale-zero Garfield timing reference failed on side "
          << side << ": "
          << PHGarfield::ForwardDriftStatusName(reference_result.status)
          << "\n";
      return false;
    }
    *reference_time[side] = reference_result.drift_time_ns;
    const double mismatch =
        std::abs(reference_result.drift_time_ns -
                 nominal_maximum_drift_time_ns) /
        nominal_maximum_drift_time_ns;
    summary.scale_zero_maximum_time_mismatch =
        std::max(
            summary.scale_zero_maximum_time_mismatch,
            mismatch);
  }
  log << "Scale-zero Garfield reference times: side0="
      << summary.scale_zero_reference_time_side0_ns
      << " ns, side1="
      << summary.scale_zero_reference_time_side1_ns
      << " ns; maximum fractional mismatch to declared PHG4 velocity="
      << summary.scale_zero_maximum_time_mismatch << "\n";
  if (summary.scale_zero_maximum_time_mismatch >
      config.maximum_nominal_time_mismatch)
  {
    error
        << "Scale-zero Garfield/PHG4 timing mismatch exceeds "
        << config.maximum_nominal_time_mismatch
        << "; make the gas, nominal E, and PHG4 calibration velocity "
           "consistent or explicitly select a larger allowed mismatch\n";
    return false;
  }

  // Capture payload content before the long integration. This identifies the
  // bytes that PHGarfield had loaded when map generation began.
  const std::string electric_field_md5 =
      PayloadMd5(garfield.GetElectricFieldMap());
  const std::string magnetic_field_md5 =
      PayloadMd5(garfield.GetResolvedMagneticFieldMap());
  const std::string gas_file_md5 =
      PayloadMd5(garfield.GetResolvedGasFile());
  if ((!garfield.GetElectricFieldMap().empty() &&
       electric_field_md5.empty()) ||
      (!garfield.GetZeroMagneticField() &&
       magnetic_field_md5.empty()) ||
      gas_file_md5.empty())
  {
    error << "Failed to fingerprint one or more PHGarfield input payloads\n";
    return false;
  }

  Provenance extended_provenance;
  extended_provenance.field_model = has_3d ? "per_side_3d" : (garfield.GetElectricFieldMap().empty() ? "nominal" : "axisymmetric_2d");
  extended_provenance.frame_field_map = garfield.GetFrameElectricFieldMap();
  extended_provenance.frame_field_map_md5 = PayloadMd5(extended_provenance.frame_field_map);
  if (!extended_provenance.frame_field_map.empty() && extended_provenance.frame_field_map_md5.empty())
  {
    error << "Failed to fingerprint the frame-field payload\n";
    return false;
  }
  for (std::size_t side = 0; side < 2; ++side)
  {
    extended_provenance.electric_field_map_3d[side] = garfield.GetElectricFieldMap3D(side);
    extended_provenance.electric_field_map_3d_md5[side] = PayloadMd5(extended_provenance.electric_field_map_3d[side]);
    extended_provenance.frame_field_map_3d[side] = garfield.GetFrameElectricFieldMap3D(side);
    extended_provenance.frame_field_map_3d_md5[side] = PayloadMd5(extended_provenance.frame_field_map_3d[side]);
    extended_provenance.frame_charge_scale[side] = garfield.GetFrameChargeScale(side);
    extended_provenance.cm_field_factor[side] = garfield.GetCMFieldRunFactor(side);
    if ((has_3d && (extended_provenance.electric_field_map_3d[side].empty() ||
                    extended_provenance.electric_field_map_3d_md5[side].empty() ||
                    !garfield.HasLoadedElectricFieldMap3D(side))) ||
        (!extended_provenance.frame_field_map_3d[side].empty() &&
         (extended_provenance.frame_field_map_3d_md5[side].empty() ||
          !garfield.HasLoadedFrameElectricFieldMap3D(side))) ||
        !std::isfinite(extended_provenance.frame_charge_scale[side]) ||
        !std::isfinite(extended_provenance.cm_field_factor[side]) || extended_provenance.cm_field_factor[side] <= 0.0)
    {
      error << "Missing, unloaded or unidentifiable side-specific field payload/settings\n";
      return false;
    }
  }
  extended_provenance.use_ifc_voltage_distortion = garfield.GetUseIFCVoltageDistortion();
  extended_provenance.use_ofc_voltage_distortion = garfield.GetUseOFCVoltageDistortion();
  extended_provenance.field_cage_offsets_v = {{garfield.GetIFCVoltageOffsetSide0(), garfield.GetIFCVoltageOffsetSide1(),
                                               garfield.GetOFCVoltageOffsetSide0(), garfield.GetOFCVoltageOffsetSide1()}};

  provenance = extended_provenance;
  provenance.nominal_electric_field_vcm = garfield.GetCMVoltageDefault();
  provenance.space_charge_scale_side0 = garfield.GetSpaceChargeScaleSide0();
  provenance.space_charge_scale_side1 = garfield.GetSpaceChargeScaleSide1();
  provenance.magnetic_field_scale = garfield.GetMagneticFieldMapScale();
  provenance.zero_magnetic_field = garfield.GetZeroMagneticField();
  provenance.electric_field_map = garfield.GetElectricFieldMap();
  provenance.magnetic_field_map = garfield.GetResolvedMagneticFieldMap();
  provenance.gas_file = garfield.GetResolvedGasFile();
  provenance.electric_field_map_md5 = electric_field_md5;
  provenance.magnetic_field_map_md5 = magnetic_field_md5;
  provenance.gas_file_md5 = gas_file_md5;
  provenance.maximum_integration_steps = PHGarfield::MaxForwardDriftSteps;
  provenance.reach_probe_scheme = config.probe_cell_centers ? "corners+center+faces+cage-biased-v1" : "corners-only-debug";
  provenance.scale_zero_reference_time_side0_ns = summary.scale_zero_reference_time_side0_ns;
  provenance.scale_zero_reference_time_side1_ns = summary.scale_zero_reference_time_side1_ns;
  provenance.scale_zero_maximum_time_mismatch = summary.scale_zero_maximum_time_mismatch;
  return true;
}

bool PHGarfieldForwardMap::WriteMetadata(TDirectory* directory,
                                         const Config& config, const Provenance& provenance, const Summary& summary,
                                         int schema_version, const std::string& schema_name)
{
  if (!directory)
  {
    return false;
  }
  TDirectory::TContext context(directory);
  const std::string canonical_text = BuildCanonicalConfiguration(config, provenance, schema_version);
  write_named("schema_name", schema_name.empty() ? (config.phi_chunk_count == 1 ? SchemaName : "PHG4TpcGarfieldFullTransportLUTShard") : schema_name);
  write_parameter("generation_phi_chunk_index", static_cast<int>(config.phi_chunk_index));
  write_parameter("generation_phi_chunk_count", static_cast<int>(config.phi_chunk_count));
  write_parameter("schema_version", schema_version);
  if (schema_version >= NumericalSchemaVersion)
  {
    write_parameter("interpolation_tolerance_um", config.interpolation_tolerance_um);
    write_named("interpolation_error_probe_scheme", schema_version >= 4 ? provenance.reach_probe_scheme : "center_and_face_centers_transverse_time_equivalent_v1");
  }
  write_named("transport_contract", TransportContract);
  write_named("coordinate_frame", "TPC_ENVELOPE_LOCAL");
  write_named("axis_order_units", "phi_rad,r_cm,z_cm");
  write_named(
      "transverse_quantity",
      "delta_r_cm,r0_delta_phi_cm");
  write_named("time_quantity", "absolute_drift_duration_ns");
  write_named("integrator", "explicit_euler");
  write_named(
      "magnetic_field_domain_policy",
      MagneticFieldDomainPolicy);
  write_named(
      "nominal_time_reference_policy",
      NominalTimeReferencePolicy);
  write_named(
      "reach_encoding",
      schema_version >= NumericalSchemaVersion ? "0=AllLost,1=Mixed,2=AllReach,3=NumericallyUnresolved" : "0=AllLost,1=Mixed,2=AllReach");
  write_named("reach_probe_scheme", provenance.reach_probe_scheme);
  write_named("configuration_label", config.configuration_label);
  // Descriptive for 2D as well, without changing its established fingerprint.
  if (provenance.field_model == "per_side_3d" || provenance.field_model == "axisymmetric_2d")
  {
    write_named("space_charge_z_coordinate_policy", SpaceChargeZCoordinatePolicy);
  }
  if (schema_version >= ExtendedSchemaVersion)
  {
    write_named("field_model", provenance.field_model);
    write_named("frame_field_map", provenance.frame_field_map);
    write_named("frame_field_map_md5", provenance.frame_field_map_md5);
    for (std::size_t side = 0; side < 2; ++side)
    {
      const std::string side_tag = "_side" + std::to_string(side);
      write_named(("electric_field_map_3d" + side_tag).c_str(), provenance.electric_field_map_3d[side]);
      write_named(("electric_field_map_3d_md5" + side_tag).c_str(), provenance.electric_field_map_3d_md5[side]);
      write_named(("frame_field_map_3d" + side_tag).c_str(), provenance.frame_field_map_3d[side]);
      write_named(("frame_field_map_3d_md5" + side_tag).c_str(), provenance.frame_field_map_3d_md5[side]);
      write_parameter(("frame_charge_scale" + side_tag).c_str(), provenance.frame_charge_scale[side]);
      write_parameter(("cm_field_factor" + side_tag).c_str(), provenance.cm_field_factor[side]);
    }
    write_parameter("use_ifc_voltage_distortion", static_cast<int>(provenance.use_ifc_voltage_distortion));
    write_parameter("use_ofc_voltage_distortion", static_cast<int>(provenance.use_ofc_voltage_distortion));
    for (std::size_t index = 0; index < 4; ++index)
    {
      write_parameter(("field_cage_offset_v_" + std::to_string(index)).c_str(), provenance.field_cage_offsets_v[index]);
    }
  }
  write_named("electric_field_map", provenance.electric_field_map);
  write_named(
      "magnetic_field_map",
      provenance.magnetic_field_map);
  write_named("gas_file", provenance.gas_file);
  write_named("electric_field_map_md5", provenance.electric_field_map_md5);
  write_named("magnetic_field_map_md5", provenance.magnetic_field_map_md5);
  write_named("gas_file_md5", provenance.gas_file_md5);
  write_named("cdb_global_tag", config.cdb_global_tag);
  write_named("canonical_configuration", canonical_text);
  write_named(
      "configuration_fingerprint_fnv1a64",
      Fingerprint(canonical_text));

  write_parameter("phi_nodes", static_cast<int>(config.phi_nodes));
  write_parameter(
      "radial_nodes",
      static_cast<int>(config.radial_nodes));
  write_parameter(
      "longitudinal_nodes",
      static_cast<int>(config.longitudinal_nodes));
  write_parameter("gas_inner_radius_cm", config.inner_radius_cm);
  write_parameter("gas_outer_radius_cm", config.outer_radius_cm);
  write_parameter(
      "central_membrane_half_width_cm",
      config.central_membrane_half_width_cm);
  write_parameter("readout_abs_z_cm", config.readout_abs_z_cm);
  write_parameter("integration_step_ns", config.integration_step_ns);
  write_parameter(
      "nominal_phg4_drift_velocity_cm_per_ns",
      config.nominal_phg4_drift_velocity_cm_per_ns);
  write_parameter(
      "maximum_nominal_time_mismatch",
      config.maximum_nominal_time_mismatch);
  write_parameter(
      "scale_zero_reference_time_side0_ns",
      summary.scale_zero_reference_time_side0_ns);
  write_parameter(
      "scale_zero_reference_time_side1_ns",
      summary.scale_zero_reference_time_side1_ns);
  write_parameter(
      "scale_zero_maximum_time_mismatch",
      summary.scale_zero_maximum_time_mismatch);
  write_parameter("tpc_place_x_cm", config.tpc_place_x_cm);
  write_parameter("tpc_place_y_cm", config.tpc_place_y_cm);
  write_parameter("tpc_place_z_cm", config.tpc_place_z_cm);
  write_parameter("tpc_rotation_x_rad", config.tpc_rotation_x_rad);
  write_parameter("tpc_rotation_y_rad", config.tpc_rotation_y_rad);
  write_parameter("tpc_rotation_z_rad", config.tpc_rotation_z_rad);
  write_parameter(
      "nominal_electric_field_vcm",
      provenance.nominal_electric_field_vcm);
  write_parameter(
      "space_charge_scale_side0",
      provenance.space_charge_scale_side0);
  write_parameter(
      "space_charge_scale_side1",
      provenance.space_charge_scale_side1);
  write_parameter(
      "magnetic_field_scale",
      provenance.magnetic_field_scale);
  write_parameter(
      "zero_magnetic_field",
      provenance.zero_magnetic_field ? 1 : 0);
  write_parameter(
      "maximum_integration_steps",
      static_cast<int>(PHGarfield::MaxForwardDriftSteps));
  write_parameter(
      "maximum_drift_time_ns",
      summary.maximum_drift_time_ns);
  write_parameter(
      "cell_center_probes",
      config.probe_cell_centers ? 1 : 0);
  write_named(
      "reach_probe_trajectories",
      std::to_string(summary.probe_trajectories));
  write_named(
      "cdb_timestamp",
      std::to_string(config.cdb_timestamp));

  return true;
}

bool PHGarfieldForwardMap::Generate(
    const PHGarfield& garfield,
    const Config& config,
    Summary& summary,
    std::ostream& log,
    std::ostream& error)
{
  Provenance provenance;
  if (!Prepare(garfield, config, provenance, summary, log, error))
  {
    return false;
  }
  const bool has_3d = !garfield.GetElectricFieldMap3D(0).empty() || !garfield.GetElectricFieldMap3D(1).empty();
  const bool has_frames = !garfield.GetFrameElectricFieldMap().empty() ||
                          !garfield.GetFrameElectricFieldMap3D(0).empty() || !garfield.GetFrameElectricFieldMap3D(1).empty();
  const int schema_version = config.interpolation_tolerance_um > 0. ? NumericalSchemaVersion : (has_3d || has_frames || garfield.GetCMFieldRunFactor(0) != 1.0 || garfield.GetCMFieldRunFactor(1) != 1.0 || garfield.GetUseIFCVoltageDistortion() || garfield.GetUseOFCVoltageDistortion() ? ExtendedSchemaVersion : SchemaVersion);
  const std::filesystem::path final_path(config.output_filename);
  std::error_code filesystem_error;
  if (!config.overwrite &&
      std::filesystem::exists(final_path, filesystem_error))
  {
    error << "Output ROOT file already exists: " << config.output_filename
          << " (use --force to replace it)\n";
    return false;
  }
  if (filesystem_error)
  {
    error << "Unable to inspect output ROOT path " << config.output_filename
          << ": " << filesystem_error.message() << "\n";
    return false;
  }

  const auto unique_stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  std::filesystem::path temporary_path;
  for (unsigned int attempt = 0; attempt < 100; ++attempt)
  {
    temporary_path =
        final_path.string() + ".partial." + std::to_string(unique_stamp) +
        "." + std::to_string(attempt);
    filesystem_error.clear();
    if (!std::filesystem::exists(temporary_path, filesystem_error) &&
        !filesystem_error)
    {
      break;
    }
    temporary_path.clear();
  }
  if (temporary_path.empty())
  {
    error << "Unable to reserve a temporary output path beside "
          << config.output_filename << "\n";
    return false;
  }

  std::unique_ptr<TFile> output(
      TFile::Open(temporary_path.string().c_str(), "CREATE"));
  if (!output || output->IsZombie())
  {
    error << "Unable to create temporary ROOT file "
          << temporary_path.string() << "\n";
    output.reset();
    std::error_code remove_error;
    std::filesystem::remove(temporary_path, remove_error);
    return false;
  }
  const auto discard_temporary_output =
      [&output, &temporary_path]()
  {
    if (output)
    {
      output->Close();
      output.reset();
    }
    std::error_code remove_error;
    std::filesystem::remove(temporary_path, remove_error);
  };

  const unsigned int phi_begin = static_cast<unsigned int>(
      static_cast<unsigned long long>(config.phi_nodes) * config.phi_chunk_index / config.phi_chunk_count);
  const unsigned int phi_end = static_cast<unsigned int>(
      static_cast<unsigned long long>(config.phi_nodes) * (config.phi_chunk_index + 1U) / config.phi_chunk_count);
  const auto owns_node = [&](const unsigned int iphi)
  {
    // One adjacent node is needed for the last owned cell, including the seam.
    return (iphi >= phi_begin && iphi < phi_end) || iphi == phi_end % config.phi_nodes;
  };
  const UniformNodeAxis phi_axis = make_phi_axis(config.phi_nodes);
  const UniformNodeAxis radial_axis = make_node_axis(
      config.radial_nodes,
      config.inner_radius_cm,
      config.outer_radius_cm);

  const std::array<UniformNodeAxis, 2> longitudinal_axes = {
      make_node_axis(
          config.longitudinal_nodes,
          -config.readout_abs_z_cm,
          -config.central_membrane_half_width_cm),
      make_node_axis(
          config.longitudinal_nodes,
          config.central_membrane_half_width_cm,
          config.readout_abs_z_cm)};

  const std::array<std::string, 2> suffix = {"negz", "posz"};
  std::array<std::unique_ptr<TH3D>, 2> delta_r;
  std::array<std::unique_ptr<TH3D>, 2> delta_rphi;
  std::array<std::unique_ptr<TH3D>, 2> drift_time;
  std::array<std::unique_ptr<TH3I>, 2> reaches_readout;
  std::array<std::unique_ptr<TH3I>, 2> reach_cell_class;
  std::array<std::unique_ptr<TH3D>, 2> interpolation_transverse_error;
  std::array<std::unique_ptr<TH3D>, 2> interpolation_time_error;

  for (unsigned int side = 0; side < 2; ++side)
  {
    const auto make_node_histogram =
        [&](const std::string& name, const std::string& title)
    {
      auto histogram = std::make_unique<TH3D>(
          name.c_str(),
          title.c_str(),
          static_cast<int>(config.phi_nodes) + 2,
          phi_axis.histogram_min,
          phi_axis.histogram_max,
          static_cast<int>(config.radial_nodes) + 2,
          radial_axis.histogram_min,
          radial_axis.histogram_max,
          static_cast<int>(config.longitudinal_nodes) + 2,
          longitudinal_axes[side].histogram_min,
          longitudinal_axes[side].histogram_max);
      histogram->SetDirectory(nullptr);
      return histogram;
    };

    delta_r[side] = make_node_histogram(
        "hDeltaR_" + suffix[side],
        "Garfield full-forward #Delta r;#phi [rad];r [cm];z [cm]");
    delta_rphi[side] = make_node_histogram(
        "hDeltaRPhi_" + suffix[side],
        "Garfield full-forward r_{0}#Delta#phi;#phi [rad];r [cm];z [cm]");
    drift_time[side] = make_node_histogram(
        "hDriftTime_" + suffix[side],
        "Garfield absolute drift duration;#phi [rad];r [cm];z [cm]");

    reaches_readout[side] = std::make_unique<TH3I>(
        ("hReachesReadoutNode_" + suffix[side]).c_str(),
        "Garfield physical node reachability;#phi [rad];r [cm];z [cm]",
        static_cast<int>(config.phi_nodes) + 2,
        phi_axis.histogram_min,
        phi_axis.histogram_max,
        static_cast<int>(config.radial_nodes) + 2,
        radial_axis.histogram_min,
        radial_axis.histogram_max,
        static_cast<int>(config.longitudinal_nodes) + 2,
        longitudinal_axes[side].histogram_min,
        longitudinal_axes[side].histogram_max);
    reaches_readout[side]->SetDirectory(nullptr);

    reach_cell_class[side] = std::make_unique<TH3I>(
        ("hReachCellClass_" + suffix[side]).c_str(),
        "Garfield reach cell class (0 lost, 1 mixed, 2 reach);"
        "#phi [rad];r [cm];z [cm]",
        static_cast<int>(config.phi_nodes),
        0.0,
        TwoPi,
        static_cast<int>(config.radial_nodes) - 1,
        config.inner_radius_cm,
        config.outer_radius_cm,
        static_cast<int>(config.longitudinal_nodes) - 1,
        longitudinal_axes[side].physical_min,
        longitudinal_axes[side].physical_max);
    reach_cell_class[side]->SetDirectory(nullptr);
    if (schema_version >= NumericalSchemaVersion)
    {
      const auto make_error = [&](const std::string& stem)
      {
        auto h = std::make_unique<TH3D>((stem + suffix[side]).c_str(), "Maximum independently probed interpolation error [um];phi;r;z",
                                        static_cast<int>(config.phi_nodes), 0., TwoPi, static_cast<int>(config.radial_nodes) - 1, config.inner_radius_cm, config.outer_radius_cm,
                                        static_cast<int>(config.longitudinal_nodes) - 1, longitudinal_axes[side].physical_min, longitudinal_axes[side].physical_max);
        h->SetDirectory(nullptr);
        return h;
      };
      interpolation_transverse_error[side] = make_error("hInterpolationTransverseError_");
      interpolation_time_error[side] = make_error("hInterpolationTimeEquivalentError_");
    }

    log << "Generating PHGarfield source nodes for "
        << (side == 0 ? "negative" : "positive") << " z\n";
    for (unsigned int iphi = 0; iphi < config.phi_nodes; ++iphi)
    {
      if (!owns_node(iphi))
      {
        continue;
      }
      const double phi_rad =
          static_cast<double>(iphi) * phi_axis.spacing;
      const int phi_bin = static_cast<int>(iphi) + 2;
      for (unsigned int ir = 0; ir < config.radial_nodes; ++ir)
      {
        const double r_cm = physical_node_value(radial_axis, ir);
        const int radial_bin = static_cast<int>(ir) + 2;
        for (unsigned int iz = 0; iz < config.longitudinal_nodes; ++iz)
        {
          const double z_cm =
              physical_node_value(longitudinal_axes[side], iz);
          const int longitudinal_bin = static_cast<int>(iz) + 2;

          NodeValue value;
          if (!run_drift(
                  garfield,
                  r_cm,
                  phi_rad,
                  z_cm,
                  config,
                  value,
                  error))
          {
            discard_temporary_output();
            return false;
          }

          ++summary.source_nodes;
          if (value.reaches_readout)
          {
            ++summary.reaching_nodes;
            summary.maximum_drift_time_ns =
                std::max(
                    summary.maximum_drift_time_ns,
                    value.drift_time_ns);
          }
          else
          {
            ++summary.lost_nodes;
          }

          delta_r[side]->SetBinContent(
              phi_bin, radial_bin, longitudinal_bin, value.delta_r_cm);
          delta_rphi[side]->SetBinContent(
              phi_bin,
              radial_bin,
              longitudinal_bin,
              value.delta_rphi_cm);
          drift_time[side]->SetBinContent(
              phi_bin,
              radial_bin,
              longitudinal_bin,
              value.drift_time_ns);
          reaches_readout[side]->SetBinContent(
              phi_bin,
              radial_bin,
              longitudinal_bin,
              value.reaches_readout ? 1 : 0);
        }
      }
    }

    fill_guard_nodes(
        *delta_r[side],
        config.phi_nodes,
        config.radial_nodes,
        config.longitudinal_nodes);
    fill_guard_nodes(
        *delta_rphi[side],
        config.phi_nodes,
        config.radial_nodes,
        config.longitudinal_nodes);
    fill_guard_nodes(
        *drift_time[side],
        config.phi_nodes,
        config.radial_nodes,
        config.longitudinal_nodes);
    fill_guard_nodes(
        *reaches_readout[side],
        config.phi_nodes,
        config.radial_nodes,
        config.longitudinal_nodes);

    log << "Classifying interpolation cells for "
        << (side == 0 ? "negative" : "positive") << " z\n";
    for (unsigned int iphi = phi_begin; iphi < phi_end; ++iphi)
    {
      const int phi_low_bin = static_cast<int>(iphi) + 2;
      const int phi_high_bin = phi_low_bin + 1;

      for (unsigned int ir = 0; ir + 1 < config.radial_nodes; ++ir)
      {
        const int radial_low_bin = static_cast<int>(ir) + 2;
        const int radial_high_bin = radial_low_bin + 1;

        for (unsigned int iz = 0;
             iz + 1 < config.longitudinal_nodes;
             ++iz)
        {
          const int longitudinal_low_bin = static_cast<int>(iz) + 2;
          const int longitudinal_high_bin = longitudinal_low_bin + 1;
          unsigned int reaching_corners = 0;
          double minimum_corner_delta_phi =
              std::numeric_limits<double>::infinity();
          double maximum_corner_delta_phi =
              -std::numeric_limits<double>::infinity();
          for (const int xbin : {phi_low_bin, phi_high_bin})
          {
            for (const int ybin : {radial_low_bin, radial_high_bin})
            {
              for (const int zbin :
                   {longitudinal_low_bin, longitudinal_high_bin})
              {
                if (reaches_readout[side]->GetBinContent(
                        xbin, ybin, zbin) == 1)
                {
                  ++reaching_corners;
                  const double source_radius_cm =
                      reaches_readout[side]->GetYaxis()->GetBinCenter(ybin);
                  const double delta_phi =
                      delta_rphi[side]->GetBinContent(
                          xbin, ybin, zbin) /
                      source_radius_cm;
                  minimum_corner_delta_phi =
                      std::min(minimum_corner_delta_phi, delta_phi);
                  maximum_corner_delta_phi =
                      std::max(maximum_corner_delta_phi, delta_phi);
                }
              }
            }
          }

          ReachCellClass classification = ReachCellClass::Mixed;
          if (reaching_corners == 0)
          {
            classification = ReachCellClass::AllLost;
          }
          else if (reaching_corners == 8)
          {
            classification = ReachCellClass::AllReach;
          }

          double maximum_transverse_error_um = 0.;
          double maximum_time_error_um = 0.;
          if (config.probe_cell_centers &&
              classification != ReachCellClass::Mixed)
          {
            const bool corners_reach =
                classification == ReachCellClass::AllReach;
            const double phi_low =
                static_cast<double>(iphi) * phi_axis.spacing;
            const double phi_high =
                static_cast<double>(iphi + 1) * phi_axis.spacing;
            const double r_low = physical_node_value(radial_axis, ir);
            const double r_high =
                physical_node_value(radial_axis, ir + 1);
            const double z_low =
                physical_node_value(longitudinal_axes[side], iz);
            const double z_high =
                physical_node_value(longitudinal_axes[side], iz + 1);

            const std::array<std::array<double, 3>, 7> probe_fractions = {{
                {{0.5, 0.5, 0.5}},
                {{0.0, 0.5, 0.5}},
                {{1.0, 0.5, 0.5}},
                {{0.5, 0.0, 0.5}},
                {{0.5, 1.0, 0.5}},
                {{0.5, 0.5, 0.0}},
                {{0.5, 0.5, 1.0}},
            }};

            bool probe_failed = false;
            const auto probe =
                [&](const std::array<double, 3>& fraction)
            {
              NodeValue probe_value;
              const bool ok = run_drift(
                  garfield,
                  r_low + fraction[1] * (r_high - r_low),
                  phi_low + fraction[0] * (phi_high - phi_low),
                  z_low + fraction[2] * (z_high - z_low),
                  config,
                  probe_value,
                  error);
              ++summary.probe_trajectories;
              if (!ok)
              {
                probe_failed = true;
                return false;
              }
              if (schema_version >= NumericalSchemaVersion && corners_reach && probe_value.reaches_readout)
              {
                const double probe_r = r_low + fraction[1] * (r_high - r_low);
                const double probe_phi = phi_low + fraction[0] * (phi_high - phi_low);
                const double probe_z = z_low + fraction[2] * (z_high - z_low);
                const double dr = delta_r[side]->Interpolate(probe_phi, probe_r, probe_z);
                const double ds = delta_rphi[side]->Interpolate(probe_phi, probe_r, probe_z);
                const double time = drift_time[side]->Interpolate(probe_phi, probe_r, probe_z);
                const double radius = probe_r + dr;
                const double angle = ds / probe_r;
                const double truth_radius = probe_r + probe_value.delta_r_cm;
                const double truth_angle = probe_value.delta_rphi_cm / probe_r;
                const double transverse_um = 1.e4 * std::hypot(radius * std::cos(angle) - truth_radius * std::cos(truth_angle), radius * std::sin(angle) - truth_radius * std::sin(truth_angle));
                const double time_um = 1.e4 * config.nominal_phg4_drift_velocity_cm_per_ns * std::abs(time - probe_value.drift_time_ns);
                if (!std::isfinite(transverse_um) || !std::isfinite(time_um))
                {
                  probe_failed = true;
                  return false;
                }
                maximum_transverse_error_um = std::max(maximum_transverse_error_um, transverse_um);
                maximum_time_error_um = std::max(maximum_time_error_um, time_um);
              }
              return probe_value.reaches_readout == corners_reach;
            };

            bool probes_agree = true;
            for (const auto& fraction : probe_fractions)
            {
              if (!probe(fraction))
              {
                probes_agree = false;
                break;
              }
            }

            // Add a probe just inside either radial cage. This specifically
            // targets a thin reach/loss boundary that corners and face centers
            // can miss.
            if (probes_agree &&
                (ir == 0 || ir + 2 == config.radial_nodes))
            {
              const double radial_fraction =
                  ir == 0 ? 0.1 : 0.9;
              probes_agree = probe(
                  {{0.5, radial_fraction, 0.5}});
            }

            if (!probes_agree)
            {
              if (probe_failed)
              {
                discard_temporary_output();
                return false;
              }
              classification = ReachCellClass::Mixed;
            }
          }

          if (classification == ReachCellClass::AllReach &&
              maximum_corner_delta_phi - minimum_corner_delta_phi >
                  std::numbers::pi)
          {
            error
                << "Wrapped delta-phi branch changes inside an all-reach "
                   "interpolation cell at side/iphi/ir/iz="
                << side << "/" << iphi << "/" << ir << "/" << iz
                << "\n";
            discard_temporary_output();
            return false;
          }

          if (schema_version >= NumericalSchemaVersion)
          {
            interpolation_transverse_error[side]->SetBinContent(static_cast<int>(iphi) + 1, static_cast<int>(ir) + 1, static_cast<int>(iz) + 1, maximum_transverse_error_um);
            interpolation_time_error[side]->SetBinContent(static_cast<int>(iphi) + 1, static_cast<int>(ir) + 1, static_cast<int>(iz) + 1, maximum_time_error_um);
            if (classification == ReachCellClass::AllReach && std::max(maximum_transverse_error_um, maximum_time_error_um) > config.interpolation_tolerance_um)
            {
              classification = ReachCellClass::NumericallyUnresolved;
            }
          }
          switch (classification)
          {
          case ReachCellClass::NumericallyUnresolved:
            ++summary.numerically_unresolved_cells;
            break;
          case ReachCellClass::AllLost:
            ++summary.all_lost_cells;
            break;
          case ReachCellClass::Mixed:
            ++summary.mixed_cells;
            break;
          case ReachCellClass::AllReach:
            ++summary.all_reach_cells;
            break;
          }

          reach_cell_class[side]->SetBinContent(
              static_cast<int>(iphi) + 1,
              static_cast<int>(ir) + 1,
              static_cast<int>(iz) + 1,
              static_cast<int>(classification));
        }
      }
    }
  }

  output->cd();
  TDirectory* transport_directory =
      output->mkdir("GarfieldFullTransport");
  if (!transport_directory)
  {
    error << "Failed to create GarfieldFullTransport directory\n";
    discard_temporary_output();
    return false;
  }
  transport_directory->cd();
  bool write_ok = true;
  for (unsigned int side = 0; side < 2; ++side)
  {
    write_ok = delta_r[side]->Write() > 0 && write_ok;
    write_ok = delta_rphi[side]->Write() > 0 && write_ok;
    write_ok = drift_time[side]->Write() > 0 && write_ok;
    write_ok = reaches_readout[side]->Write() > 0 && write_ok;
    write_ok = reach_cell_class[side]->Write() > 0 && write_ok;
    if (schema_version >= NumericalSchemaVersion)
    {
      write_ok = interpolation_transverse_error[side]->Write() > 0 && write_ok;
      write_ok = interpolation_time_error[side]->Write() > 0 && write_ok;
    }
  }

  TDirectory* metadata_directory =
      transport_directory->mkdir("Metadata");
  if (!metadata_directory)
  {
    error << "Failed to create GarfieldFullTransport/Metadata directory\n";
    discard_temporary_output();
    return false;
  }
  metadata_directory->cd();

  if (!WriteMetadata(metadata_directory, config, provenance, summary, schema_version))
  {
    discard_temporary_output();
    return false;
  }

  output->Write();
  output->Flush();
  write_ok =
      write_ok &&
      !output->TestBit(TFile::kWriteError) &&
      output->GetBytesWritten() > 0;
  if (!write_ok)
  {
    error << "Failed while writing temporary ROOT file "
          << temporary_path.string() << "\n";
    discard_temporary_output();
    return false;
  }
  output->Close();
  output.reset();

  filesystem_error.clear();
  if (!config.overwrite &&
      std::filesystem::exists(final_path, filesystem_error))
  {
    error << "Output ROOT file appeared during generation: "
          << config.output_filename << "\n";
    std::error_code remove_error;
    std::filesystem::remove(temporary_path, remove_error);
    return false;
  }
  if (filesystem_error)
  {
    error << "Unable to recheck output path " << config.output_filename
          << ": " << filesystem_error.message() << "\n";
    std::error_code remove_error;
    std::filesystem::remove(temporary_path, remove_error);
    return false;
  }

  std::filesystem::rename(temporary_path, final_path, filesystem_error);
  if (filesystem_error)
  {
    error << "Unable to atomically publish " << config.output_filename
          << ": " << filesystem_error.message() << "\n";
    std::error_code remove_error;
    std::filesystem::remove(temporary_path, remove_error);
    return false;
  }

  log << "Wrote " << config.output_filename << " with "
      << summary.source_nodes << " source nodes ("
      << summary.reaching_nodes << " reaching, "
      << summary.lost_nodes << " lost) and "
      << summary.all_reach_cells << " all-reach, "
      << summary.all_lost_cells << " all-lost, "
      << summary.mixed_cells << " mixed cells and " << summary.numerically_unresolved_cells << " numerically unresolved cells\n";
  return true;
}

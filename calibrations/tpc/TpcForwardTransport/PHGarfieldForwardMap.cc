#include "PHGarfieldForwardMap.h"
#include "PHGarfieldFingerprint.h"

#include <TMD5.h>

#include <algorithm>
#include <array>
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

  std::string regular_file_md5(const std::filesystem::path& filename)
  {
    std::unique_ptr<TMD5> checksum(
        TMD5::FileChecksum(filename.string().c_str()));
    return checksum ? checksum->AsString() : std::string{};
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

    // Match PHGarfieldForwardTransport::InitializeGas: only ordered MERGED_E###.gas payloads
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
    const PHGarfieldForwardTransport::ForwardDriftResult& result,
    NodeValue& value)
{
  value = NodeValue{};
  if (!result.reached_readout ||
      result.status != PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout ||
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
  require(config.phi_nodes >= 3, "At least three phi nodes are required");
  require(config.radial_nodes >= 2, "At least two radial nodes are required");
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
  require(std::isfinite(config.max_radial_field_extension_cm) &&
              config.max_radial_field_extension_cm >= 0.,
          "Maximum radial field extension must be finite and nonnegative");
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
      "Reference drift-time warning threshold must be finite and nonnegative");
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
    const Provenance& provenance)
{
  std::ostringstream canonical;
  const auto append_string =
      [&canonical](const char* key, const std::string& value)
  {
    canonical << ";" << key << "[" << value.size() << "]=" << value;
  };
  canonical << std::setprecision(17)
            << "schema=" << SchemaName
            << ";version=" << SchemaVersion
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
            << 1  // Full 27-point stencil is mandatory.
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
  canonical << ";interpolation_tolerance_um=" << config.interpolation_tolerance_um;
  if (config.max_radial_field_extension_cm > 0.)
  {
    canonical << ";max_radial_field_extension_cm=" << config.max_radial_field_extension_cm;
  }
  append_string("interpolation_error_probe_scheme", provenance.reach_probe_scheme);
  return canonical.str();
}

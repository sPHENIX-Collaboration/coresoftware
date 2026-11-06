#include "PHGarfieldForwardMap.h"
#include "PHGarfieldFingerprint.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TMD5.h>
#include <TNamed.h>
#include <TParameter.h>

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
            << 1  // Frozen full-stencil identity, shared with existing adaptive tables.
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
  provenance.reach_probe_scheme = "corners+center+faces+cage-biased-v1";
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
  write_named("schema_name", schema_name.empty() ? SchemaName : schema_name);
  write_parameter("generation_phi_chunk_index", 0);
  write_parameter("generation_phi_chunk_count", 1);
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
      1);
  write_named(
      "reach_probe_trajectories",
      std::to_string(summary.probe_trajectories));
  write_named(
      "cdb_timestamp",
      std::to_string(config.cdb_timestamp));

  return true;
}

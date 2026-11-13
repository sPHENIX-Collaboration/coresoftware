#include "../PHGarfieldForwardMap.h"

#include <TMemFile.h>
#include <TNamed.h>
#include <TParameter.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <utility>

namespace
{
  bool close(const double first, const double second)
  {
    return std::abs(first - second) < 1.0e-12;
  }

  bool expect_config(
      const PHGarfieldForwardMap::Config& config,
      const std::string& description, bool valid)
  {
    std::ostringstream error;
    if (PHGarfieldForwardMap::ValidateConfig(config, error) == valid &&
        (valid || !error.str().empty()))
    {
      return true;
    }

    std::cerr << "Unexpected ValidateConfig result for " << description << ": " << error.str() << '\n';
    return false;
  }
}  // namespace

int main()
{
  PHGarfieldForwardMap::Config valid_config;
  valid_config.output_filename = "unused-forward-map.root";
  if (!expect_config(valid_config, "the default configuration", true))
  {
    return 1;
  }

  using Config = PHGarfieldForwardMap::Config;
  const double infinity = std::numeric_limits<double>::infinity();
  const std::pair<const char*, void (*)(Config&)> rejected_configs[] = {
      {"an empty output filename", [](auto& config)
       { config.output_filename.clear(); }},
      {"too few phi nodes", [](auto& config)
       { config.phi_nodes = 2; }},
      {"too few radial nodes", [](auto& config)
       { config.radial_nodes = 1; }},
      {"too few longitudinal nodes", [](auto& config)
       { config.longitudinal_nodes = 1; }},
      {"a wrapped-like phi node count", [](auto& config)
       { config.phi_nodes = std::numeric_limits<unsigned int>::max(); }},
      {"a wrapped-like radial node count", [](auto& config)
       { config.radial_nodes = std::numeric_limits<unsigned int>::max(); }},
      {"a wrapped-like longitudinal node count", [](auto& config)
       { config.longitudinal_nodes = std::numeric_limits<unsigned int>::max(); }},
      {"a negative interpolation tolerance", [](auto& config)
       { config.interpolation_tolerance_um = -1.; }},
      {"a nonfinite interpolation tolerance", [](auto& config)
       { config.interpolation_tolerance_um = std::numeric_limits<double>::quiet_NaN(); }},
      {"an excessive total source-node count", [](auto& config)
       { config.phi_nodes = 1000; config.radial_nodes = 1000; config.longitudinal_nodes = 6; }},
      {"a zero inner radius", [](auto& config)
       { config.inner_radius_cm = 0.0; }},
      {"a nonfinite inner radius", [](auto& config)
       { config.inner_radius_cm = std::numeric_limits<double>::quiet_NaN(); }},
      {"an outer radius equal to the inner radius", [](auto& config)
       { config.outer_radius_cm = config.inner_radius_cm; }},
      {"a nonfinite outer radius", [](auto& config)
       { config.outer_radius_cm = std::numeric_limits<double>::infinity(); }},
      {"a negative radial field extension", [](auto& config)
       { config.max_radial_field_extension_cm = -1.; }},
      {"a nonfinite radial field extension", [](auto& config)
       { config.max_radial_field_extension_cm = std::numeric_limits<double>::quiet_NaN(); }},
      {"a zero central-membrane half-width", [](auto& config)
       { config.central_membrane_half_width_cm = 0.0; }},
      {"a nonfinite central-membrane half-width", [](auto& config)
       { config.central_membrane_half_width_cm = std::numeric_limits<double>::quiet_NaN(); }},
      {"a readout plane at the central membrane", [](auto& config)
       { config.readout_abs_z_cm = config.central_membrane_half_width_cm; }},
      {"a nonfinite readout plane", [](auto& config)
       { config.readout_abs_z_cm = std::numeric_limits<double>::infinity(); }},
      {"a nonfinite TPC placement", [](auto& config)
       { config.tpc_place_y_cm = std::numeric_limits<double>::quiet_NaN(); }},
      {"a nonfinite TPC rotation", [](auto& config)
       { config.tpc_rotation_z_rad = std::numeric_limits<double>::infinity(); }},
      {"a zero integration step", [](auto& config)
       { config.integration_step_ns = 0.0; }},
      {"a negative integration step", [](auto& config)
       { config.integration_step_ns = -1.0; }},
      {"a nonfinite integration step", [](auto& config)
       { config.integration_step_ns = std::numeric_limits<double>::quiet_NaN(); }},
      {"a zero nominal PHG4 drift velocity", [](auto& config)
       { config.nominal_phg4_drift_velocity_cm_per_ns = 0.0; }},
      {"a negative nominal PHG4 drift velocity", [](auto& config)
       { config.nominal_phg4_drift_velocity_cm_per_ns = -0.001; }},
      {"a nonfinite nominal PHG4 drift velocity", [](auto& config)
       { config.nominal_phg4_drift_velocity_cm_per_ns = std::numeric_limits<double>::infinity(); }},
      {"a negative nominal-time mismatch", [](auto& config)
       { config.maximum_nominal_time_mismatch = -0.001; }},
      {"a nonfinite nominal-time mismatch", [](auto& config)
       { config.maximum_nominal_time_mismatch = std::numeric_limits<double>::quiet_NaN(); }},
  };
  for (const auto& [description, mutate] : rejected_configs)
  {
    auto config = valid_config;
    mutate(config);
    if (!expect_config(config, description, false))
    {
      return 1;
    }
  }

  constexpr double r0 = 40.0;
  constexpr double phi0 = 2.0 * std::numbers::pi - 0.01;
  constexpr double endpoint_radius = 41.25;
  constexpr double endpoint_phi = 0.02;

  PHGarfieldForwardTransport::ForwardDriftResult result;
  result.status = PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout;
  result.reached_readout = true;
  result.x_cm = endpoint_radius * std::cos(endpoint_phi);
  result.y_cm = endpoint_radius * std::sin(endpoint_phi);
  result.drift_time_ns = 1234.5;

  PHGarfieldForwardMap::NodeValue encoded;
  if (!PHGarfieldForwardMap::EncodeNode(r0, phi0, result, encoded))
  {
    std::cerr << "EncodeNode rejected a valid endpoint\n";
    return 1;
  }

  const double reconstructed_radius = r0 + encoded.delta_r_cm;
  const double reconstructed_phi = phi0 + encoded.delta_rphi_cm / r0;
  if (!encoded.reaches_readout ||
      !close(reconstructed_radius, endpoint_radius) ||
      !close(
          std::remainder(
              reconstructed_phi - endpoint_phi,
              2.0 * std::numbers::pi),
          0.0) ||
      !close(encoded.drift_time_ns, result.drift_time_ns))
  {
    std::cerr << "Full-D cylindrical encoding did not round-trip\n";
    return 1;
  }

  result.reached_readout = false;
  if (PHGarfieldForwardMap::EncodeNode(r0, phi0, result, encoded))
  {
    std::cerr << "EncodeNode accepted a lost endpoint\n";
    return 1;
  }

  PHGarfieldForwardTransport::ForwardDriftResult readout_result;
  readout_result.status =
      PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout;
  readout_result.reached_readout = true;
  readout_result.x_cm = r0 * std::cos(phi0);
  readout_result.y_cm = r0 * std::sin(phi0);
  readout_result.drift_time_ns = 0.0;
  if (!PHGarfieldForwardMap::EncodeNode(
          r0, phi0, readout_result, encoded) ||
      !encoded.reaches_readout ||
      !close(encoded.delta_r_cm, 0.0) ||
      !close(encoded.delta_rphi_cm, 0.0) ||
      !close(encoded.drift_time_ns, 0.0))
  {
    std::cerr << "EncodeNode rejected an exact readout endpoint\n";
    return 1;
  }

  using Result = PHGarfieldForwardTransport::ForwardDriftResult;
  const std::pair<const char*, void (*)(Result&)> rejected_endpoints[] = {
      {"an inconsistent reaching status", [](auto& value)
       { value.status = PHGarfieldForwardTransport::ForwardDriftStatus::LeftDriftVolume; }},
      {"a nonfinite endpoint", [](auto& value)
       { value.x_cm = std::numeric_limits<double>::quiet_NaN(); }},
      {"a negative drift time", [](auto& value)
       { value.drift_time_ns = -1.0; }},
  };
  for (const auto& [description, mutate] : rejected_endpoints)
  {
    auto invalid_result = readout_result;
    mutate(invalid_result);
    if (PHGarfieldForwardMap::EncodeNode(r0, phi0, invalid_result, encoded))
    {
      std::cerr << "EncodeNode accepted " << description << '\n';
      return 1;
    }
  }

  if (PHGarfieldForwardMap::EncodeNode(
          0.0, phi0, readout_result, encoded) ||
      PHGarfieldForwardMap::EncodeNode(
          r0, infinity, readout_result, encoded))
  {
    std::cerr << "EncodeNode accepted an invalid source coordinate\n";
    return 1;
  }

  PHGarfieldForwardMap::Provenance first_provenance;
  first_provenance.electric_field_map =
      "foo;magnetic_field_map=bar";
  first_provenance.magnetic_field_map = "baz";
  auto second_provenance = first_provenance;
  second_provenance.electric_field_map = "foo";
  second_provenance.magnetic_field_map =
      "bar;magnetic_field_map=baz";
  if (PHGarfieldForwardMap::BuildCanonicalConfiguration(
          valid_config, first_provenance) ==
      PHGarfieldForwardMap::BuildCanonicalConfiguration(
          valid_config, second_provenance))
  {
    std::cerr
        << "Canonical metadata encoding is ambiguous for delimiter text\n";
    return 1;
  }

  PHGarfieldForwardMap::Provenance extended;
  extended.field_model = "per_side_3d";
  extended.electric_field_map_3d = {{"south.root", "north.root"}};
  extended.electric_field_map_3d_md5 = {{std::string(32, 'a'), std::string(32, 'b')}};
  extended.frame_field_map_3d = {{"frames_south.root", "frames_north.root"}};
  extended.frame_field_map_3d_md5 = {{std::string(32, 'c'), std::string(32, 'd')}};
  extended.frame_charge_scale = {{-180.0, -180.0}};
  const auto canonical = [&](const auto& provenance)
  {
    return PHGarfieldForwardMap::BuildCanonicalConfiguration(valid_config, provenance);
  };
  const auto reference = canonical(extended);
  auto strict_config = valid_config;
  strict_config.max_radial_field_extension_cm = 0.;
  const auto strict_canonical = PHGarfieldForwardMap::BuildCanonicalConfiguration(strict_config, extended);
  if (valid_config.max_radial_field_extension_cm != .3 ||
      reference == strict_canonical ||
      reference.find("max_radial_field_extension_cm=") == std::string::npos ||
      strict_canonical.find("max_radial_field_extension_cm=") != std::string::npos)
  {
    std::cerr << "Radial field extension is not bound to configuration identity\n";
    return 1;
  }
  if (reference.find(PHGarfieldForwardMap::SpaceChargeZCoordinatePolicy) == std::string::npos ||
      canonical(first_provenance).find("space_charge_z_coordinate_policy") != std::string::npos)
  {
    std::cerr << "3D coordinate policy is unbound\n";
    return 1;
  }

  const auto check_bound = [&](const auto& mutate)
  {
    auto changed = extended;
    mutate(changed);
    return canonical(changed) != reference;
  };
  if (!check_bound([](auto& p)
                   { p.electric_field_map_3d_md5[0][0] = 'e'; }) ||
      !check_bound([](auto& p)
                   { p.electric_field_map_3d_md5[1][0] = 'e'; }) ||
      !check_bound([](auto& p)
                   { p.frame_field_map_3d_md5[0][0] = 'e'; }) ||
      !check_bound([](auto& p)
                   { p.frame_field_map_3d_md5[1][0] = 'e'; }) ||
      !check_bound([](auto& p)
                   { p.frame_charge_scale[0] = -179.0; }) ||
      !check_bound([](auto& p)
                   { p.frame_charge_scale[1] = -179.0; }) ||
      !check_bound([](auto& p)
                   { p.cm_field_factor[1] = 1.01; }) ||
      !check_bound([](auto& p)
                   { p.field_cage_offsets_v[0] = 1.0; }) ||
      !check_bound([](auto& p)
                   { p.use_ifc_voltage_distortion = true; }))
  {
    std::cerr << "Physical field provenance is not bound\n";
    return 1;
  }

  PHGarfieldForwardMap::Provenance axisymmetric;
  axisymmetric.field_model = "axisymmetric_2d";
  axisymmetric.electric_field_map = "axisymmetric.root";
  TMemFile output("forward-policy.root", "RECREATE");
  auto* metadata = output.mkdir("Metadata");
  if (!PHGarfieldForwardMap::WriteMetadata(metadata, valid_config, axisymmetric, {}))
  {
    return 1;
  }
  const auto* policy = dynamic_cast<const TNamed*>(metadata->Get("space_charge_z_coordinate_policy"));
  const auto* extension = dynamic_cast<const TParameter<double>*>(metadata->Get("max_radial_field_extension_cm"));
  if (!policy || std::string(policy->GetTitle()) != PHGarfieldForwardMap::SpaceChargeZCoordinatePolicy ||
      !extension || extension->GetVal() != .3)
  {
    std::cerr << "2D coordinate metadata is missing\n";
    return 1;
  }
  return 0;
}

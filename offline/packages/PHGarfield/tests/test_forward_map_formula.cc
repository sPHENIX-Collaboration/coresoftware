#include "../PHGarfieldForwardMap.h"

#include <TMemFile.h>
#include <TNamed.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>

namespace
{
  bool close(const double first, const double second)
  {
    return std::abs(first - second) < 1.0e-12;
  }

  bool expect_valid(
      const PHGarfieldForwardMap::Config& config,
      const std::string& description)
  {
    std::ostringstream error;
    if (PHGarfieldForwardMap::ValidateConfig(config, error))
    {
      return true;
    }

    std::cerr << "ValidateConfig rejected " << description << ": "
              << error.str();
    return false;
  }

  bool expect_invalid(
      const PHGarfieldForwardMap::Config& config,
      const std::string& description)
  {
    std::ostringstream error;
    if (!PHGarfieldForwardMap::ValidateConfig(config, error) &&
        !error.str().empty())
    {
      return true;
    }

    std::cerr << "ValidateConfig accepted " << description << "\n";
    return false;
  }
}  // namespace

int main()
{
  PHGarfieldForwardMap::Config valid_config;
  valid_config.output_filename = "unused-forward-map.root";
  if (!expect_valid(valid_config, "the default configuration"))
  {
    return 1;
  }

  const auto rejects =
      [&valid_config](const std::string& description, const auto& mutate)
  {
    auto config = valid_config;
    mutate(config);
    return expect_invalid(config, description);
  };

  if (!rejects(
          "an empty output filename",
          [](auto& config)
          { config.output_filename.clear(); }) ||
      !rejects(
          "too few phi nodes",
          [](auto& config)
          { config.phi_nodes = 2; }) ||
      !rejects(
          "too few radial nodes",
          [](auto& config)
          { config.radial_nodes = 1; }) ||
      !rejects(
          "too few longitudinal nodes",
          [](auto& config)
          { config.longitudinal_nodes = 1; }) ||
      !rejects(
          "a wrapped-like phi node count",
          [](auto& config)
          {
            config.phi_nodes =
                std::numeric_limits<unsigned int>::max();
          }) ||
      !rejects(
          "a wrapped-like radial node count",
          [](auto& config)
          {
            config.radial_nodes =
                std::numeric_limits<unsigned int>::max();
          }) ||
      !rejects(
          "a wrapped-like longitudinal node count",
          [](auto& config)
          {
            config.longitudinal_nodes =
                std::numeric_limits<unsigned int>::max();
          }) ||
      !rejects("a negative interpolation tolerance", [](auto& config)
               { config.interpolation_tolerance_um = -1.; }) ||
      !rejects("a nonfinite interpolation tolerance", [](auto& config)
               { config.interpolation_tolerance_um = std::numeric_limits<double>::quiet_NaN(); }) ||
      !rejects("an unprobed numerical mask", [](auto& config)
               { config.interpolation_tolerance_um = 3.; config.probe_cell_centers = false; }) ||
      !rejects("a zero phi chunk count", [](auto& config)
               { config.phi_chunk_count = 0; }) ||
      !rejects("too many phi chunks", [](auto& config)
               { config.phi_chunk_count = config.phi_nodes + 1; }) ||
      !rejects("a phi chunk index outside the partition", [](auto& config)
               { config.phi_chunk_index = 1; }) ||
      !rejects(
          "an excessive total source-node count",
          [](auto& config)
          {
            config.phi_nodes = 1000;
            config.radial_nodes = 1000;
            config.longitudinal_nodes = 6;
          }))
  {
    return 1;
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  if (!rejects(
          "a zero inner radius",
          [](auto& config)
          { config.inner_radius_cm = 0.0; }) ||
      !rejects(
          "a nonfinite inner radius",
          [nan](auto& config)
          { config.inner_radius_cm = nan; }) ||
      !rejects(
          "an outer radius equal to the inner radius",
          [](auto& config)
          {
            config.outer_radius_cm = config.inner_radius_cm;
          }) ||
      !rejects(
          "a nonfinite outer radius",
          [infinity](auto& config)
          {
            config.outer_radius_cm = infinity;
          }) ||
      !rejects(
          "a zero central-membrane half-width",
          [](auto& config)
          {
            config.central_membrane_half_width_cm = 0.0;
          }) ||
      !rejects(
          "a nonfinite central-membrane half-width",
          [nan](auto& config)
          {
            config.central_membrane_half_width_cm = nan;
          }) ||
      !rejects(
          "a readout plane at the central membrane",
          [](auto& config)
          {
            config.readout_abs_z_cm =
                config.central_membrane_half_width_cm;
          }) ||
      !rejects(
          "a nonfinite readout plane",
          [infinity](auto& config)
          {
            config.readout_abs_z_cm = infinity;
          }) ||
      !rejects(
          "a nonfinite TPC placement",
          [nan](auto& config)
          { config.tpc_place_y_cm = nan; }) ||
      !rejects(
          "a nonfinite TPC rotation",
          [infinity](auto& config)
          {
            config.tpc_rotation_z_rad = infinity;
          }))
  {
    return 1;
  }

  if (!rejects(
          "a zero integration step",
          [](auto& config)
          { config.integration_step_ns = 0.0; }) ||
      !rejects(
          "a negative integration step",
          [](auto& config)
          { config.integration_step_ns = -1.0; }) ||
      !rejects(
          "a nonfinite integration step",
          [nan](auto& config)
          { config.integration_step_ns = nan; }) ||
      !rejects(
          "a zero nominal PHG4 drift velocity",
          [](auto& config)
          {
            config.nominal_phg4_drift_velocity_cm_per_ns = 0.0;
          }) ||
      !rejects(
          "a negative nominal PHG4 drift velocity",
          [](auto& config)
          {
            config.nominal_phg4_drift_velocity_cm_per_ns = -0.001;
          }) ||
      !rejects(
          "a nonfinite nominal PHG4 drift velocity",
          [infinity](auto& config)
          {
            config.nominal_phg4_drift_velocity_cm_per_ns = infinity;
          }) ||
      !rejects(
          "a negative nominal-time mismatch",
          [](auto& config)
          {
            config.maximum_nominal_time_mismatch = -0.001;
          }) ||
      !rejects(
          "a nonfinite nominal-time mismatch",
          [nan](auto& config)
          {
            config.maximum_nominal_time_mismatch = nan;
          }))
  {
    return 1;
  }

  constexpr double r0 = 40.0;
  constexpr double phi0 = 2.0 * std::numbers::pi - 0.01;
  constexpr double endpoint_radius = 41.25;
  constexpr double endpoint_phi = 0.02;

  PHGarfield::ForwardDriftResult result;
  result.status = PHGarfield::ForwardDriftStatus::ReachedReadout;
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

  PHGarfield::ForwardDriftResult readout_result;
  readout_result.status =
      PHGarfield::ForwardDriftStatus::ReachedReadout;
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

  auto invalid_result = readout_result;
  invalid_result.status =
      PHGarfield::ForwardDriftStatus::LeftDriftVolume;
  if (PHGarfieldForwardMap::EncodeNode(
          r0, phi0, invalid_result, encoded))
  {
    std::cerr << "EncodeNode accepted an inconsistent reaching status\n";
    return 1;
  }

  invalid_result = readout_result;
  invalid_result.x_cm = nan;
  if (PHGarfieldForwardMap::EncodeNode(
          r0, phi0, invalid_result, encoded))
  {
    std::cerr << "EncodeNode accepted a nonfinite endpoint\n";
    return 1;
  }

  invalid_result = readout_result;
  invalid_result.drift_time_ns = -1.0;
  if (PHGarfieldForwardMap::EncodeNode(
          r0, phi0, invalid_result, encoded))
  {
    std::cerr << "EncodeNode accepted a negative drift time\n";
    return 1;
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
  const auto canonical = [&](const auto& provenance, int schema)
  {
    return PHGarfieldForwardMap::BuildCanonicalConfiguration(valid_config, provenance, schema);
  };
  const auto reference = canonical(extended, PHGarfieldForwardMap::ExtendedSchemaVersion);
  if (reference.find(PHGarfieldForwardMap::SpaceChargeZCoordinatePolicy) == std::string::npos ||
      canonical(first_provenance, PHGarfieldForwardMap::SchemaVersion).find("space_charge_z_coordinate_policy") != std::string::npos)
  {
    std::cerr << "3D coordinate policy is unbound or legacy 2D identity changed\n";
    return 1;
  }

  const auto check_bound = [&](const auto& mutate)
  {
    auto changed = extended;
    mutate(changed);
    return canonical(changed, PHGarfieldForwardMap::ExtendedSchemaVersion) != reference &&
           canonical(changed, PHGarfieldForwardMap::SchemaVersion) == canonical(extended, PHGarfieldForwardMap::SchemaVersion);
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
    std::cerr << "Extended field provenance is not bound, or legacy canonical text changed\n";
    return 1;
  }

  PHGarfieldForwardMap::Provenance axisymmetric;
  axisymmetric.field_model = "axisymmetric_2d";
  axisymmetric.electric_field_map = "axisymmetric.root";
  for (const int schema : {PHGarfieldForwardMap::SchemaVersion, PHGarfieldForwardMap::ExtendedSchemaVersion})
  {
    TMemFile output("forward-policy.root", "RECREATE");
    auto* metadata = output.mkdir("Metadata");
    if (!PHGarfieldForwardMap::WriteMetadata(metadata, valid_config, axisymmetric, {}, schema))
    {
      return 1;
    }
    const auto* policy = dynamic_cast<const TNamed*>(metadata->Get("space_charge_z_coordinate_policy"));
    if (!policy || std::string(policy->GetTitle()) != PHGarfieldForwardMap::SpaceChargeZCoordinatePolicy ||
        canonical(axisymmetric, schema).find("space_charge_z_coordinate_policy") != std::string::npos)
    {
      std::cerr << "2D coordinate metadata is missing or changed its established fingerprint\n";
      return 1;
    }
  }
  TMemFile invalid_map("invalid-regular.root", "RECREATE");
  auto* metadata = invalid_map.mkdir("GarfieldFullTransport")->mkdir("Metadata");
  auto invalid_config = valid_config;
  invalid_config.integration_step_ns = -1.;
  if (!PHGarfieldForwardMap::WriteMetadata(metadata, invalid_config, {}, {}, PHGarfieldForwardMap::SchemaVersion))
  {
    return 1;
  }
  std::ostringstream errors, unexpected;
  auto* previous = std::cerr.rdbuf(unexpected.rdbuf());
  const bool accepted = PHGarfieldForwardMap::ValidateRegular(invalid_map, errors);
  std::cerr.rdbuf(previous);
  if (accepted || !unexpected.str().empty() || errors.str().find("Euler integration step") == std::string::npos ||
      errors.str().find("invalid merged geometry/numerical configuration") == std::string::npos)
  {
    std::cerr << "Regular-map validation did not preserve the caller's diagnostic stream\n";
    return 1;
  }
  return 0;
}

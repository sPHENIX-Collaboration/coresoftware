#include "PHGarfieldAdaptiveMap.h"
#include "PHGarfieldForwardMap.h"

#include <TDirectory.h>
#include <TNamed.h>
#include <TParameter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <string>
#include <utility>

namespace
{
  template <class T>
  void write_parameters(std::initializer_list<std::pair<const char*, T>> parameters)
  {
    for (const auto& [name, value] : parameters)
    {
      TParameter<T> parameter(name, value);
      parameter.Write();
    }
  }

  void write_named(std::initializer_list<std::pair<const char*, std::string>> values)
  {
    for (const auto& [name, value] : values)
    {
      TNamed named(name, value.c_str());
      named.Write();
    }
  }

}  // namespace

bool PHGarfieldForwardMap::Prepare(const PHGarfieldForwardTransport& garfield,
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
    error << "Adaptive generation requires the explicit per-side 3D model\n";
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
    PHGarfieldForwardTransport::ForwardDriftResult reference_result;
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
          << PHGarfieldForwardTransport::ForwardDriftStatusName(reference_result.status)
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
  if (has_3d)
  {
    extended_provenance.field_model = "per_side_3d";
  }
  else if (garfield.GetElectricFieldMap().empty())
  {
    extended_provenance.field_model = "nominal";
  }
  else
  {
    extended_provenance.field_model = "axisymmetric_2d";
  }
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
  provenance.maximum_integration_steps = PHGarfieldForwardTransport::MaxForwardDriftSteps;
  provenance.reach_probe_scheme = "corners+center+faces+cage-biased-v1";
  provenance.scale_zero_reference_time_side0_ns = summary.scale_zero_reference_time_side0_ns;
  provenance.scale_zero_reference_time_side1_ns = summary.scale_zero_reference_time_side1_ns;
  provenance.scale_zero_maximum_time_mismatch = summary.scale_zero_maximum_time_mismatch;
  return true;
}

bool PHGarfieldForwardMap::WriteMetadata(TDirectory* directory,
                                         const Config& config, const Provenance& provenance, const Summary& summary)
{
  if (!directory)
  {
    return false;
  }
  TDirectory::TContext context(directory);
  const std::string canonical_text = BuildCanonicalConfiguration(config, provenance);
  write_named({{"schema_name", SchemaName}});
  write_parameters<int>({{"schema_version", SchemaVersion}});
  write_parameters<double>({{"interpolation_tolerance_um", config.interpolation_tolerance_um}});
  write_named({{"interpolation_error_probe_scheme", provenance.reach_probe_scheme}});
  write_named({{"transport_contract", TransportContract},
               {"coordinate_frame", "TPC_ENVELOPE_LOCAL"},
               {"axis_order_units", "phi_rad,r_cm,z_cm"},
               {"transverse_quantity", "delta_r_cm,r0_delta_phi_cm"},
               {"time_quantity", "absolute_drift_duration_ns"},
               {"integrator", "explicit_euler"},
               {"magnetic_field_domain_policy", MagneticFieldDomainPolicy},
               {"nominal_time_reference_policy", NominalTimeReferencePolicy},
               {"reach_encoding", PHGarfieldAdaptiveMap::LeafEncoding},
               {"reach_probe_scheme", provenance.reach_probe_scheme},
               {"configuration_label", config.configuration_label}});
  // The same longitudinal coordinate convention applies to 2D and 3D maps.
  if (provenance.field_model == "per_side_3d" || provenance.field_model == "axisymmetric_2d")
  {
    write_named({{"space_charge_z_coordinate_policy", SpaceChargeZCoordinatePolicy}});
  }
  write_named({{"field_model", provenance.field_model},
               {"frame_field_map", provenance.frame_field_map},
               {"frame_field_map_md5", provenance.frame_field_map_md5}});
  for (std::size_t side = 0; side < 2; ++side)
  {
    const std::string side_tag = "_side" + std::to_string(side);
    write_named({{("electric_field_map_3d" + side_tag).c_str(), provenance.electric_field_map_3d[side]},
                 {("electric_field_map_3d_md5" + side_tag).c_str(), provenance.electric_field_map_3d_md5[side]},
                 {("frame_field_map_3d" + side_tag).c_str(), provenance.frame_field_map_3d[side]},
                 {("frame_field_map_3d_md5" + side_tag).c_str(), provenance.frame_field_map_3d_md5[side]}});
    write_parameters<double>({{("frame_charge_scale" + side_tag).c_str(), provenance.frame_charge_scale[side]},
                              {("cm_field_factor" + side_tag).c_str(), provenance.cm_field_factor[side]}});
  }
  write_parameters<int>({{"use_ifc_voltage_distortion", static_cast<int>(provenance.use_ifc_voltage_distortion)},
                         {"use_ofc_voltage_distortion", static_cast<int>(provenance.use_ofc_voltage_distortion)}});
  for (std::size_t index = 0; index < 4; ++index)
  {
    write_parameters<double>({{("field_cage_offset_v_" + std::to_string(index)).c_str(), provenance.field_cage_offsets_v[index]}});
  }
  write_named({{"electric_field_map", provenance.electric_field_map},
               {"magnetic_field_map", provenance.magnetic_field_map},
               {"gas_file", provenance.gas_file},
               {"electric_field_map_md5", provenance.electric_field_map_md5},
               {"magnetic_field_map_md5", provenance.magnetic_field_map_md5},
               {"gas_file_md5", provenance.gas_file_md5},
               {"cdb_global_tag", config.cdb_global_tag},
               {"canonical_configuration", canonical_text},
               {"configuration_fingerprint_fnv1a64", Fingerprint(canonical_text)}});

  write_parameters<int>({{"phi_nodes", static_cast<int>(config.phi_nodes)},
                         {"radial_nodes", static_cast<int>(config.radial_nodes)},
                         {"longitudinal_nodes", static_cast<int>(config.longitudinal_nodes)}});
  write_parameters<double>({{"gas_inner_radius_cm", config.inner_radius_cm},
                            {"gas_outer_radius_cm", config.outer_radius_cm},
                            {"central_membrane_half_width_cm", config.central_membrane_half_width_cm},
                            {"readout_abs_z_cm", config.readout_abs_z_cm},
                            {"integration_step_ns", config.integration_step_ns},
                            {"nominal_phg4_drift_velocity_cm_per_ns", config.nominal_phg4_drift_velocity_cm_per_ns},
                            {"maximum_nominal_time_mismatch", config.maximum_nominal_time_mismatch},
                            {"scale_zero_reference_time_side0_ns", summary.scale_zero_reference_time_side0_ns},
                            {"scale_zero_reference_time_side1_ns", summary.scale_zero_reference_time_side1_ns},
                            {"scale_zero_maximum_time_mismatch", summary.scale_zero_maximum_time_mismatch},
                            {"tpc_place_x_cm", config.tpc_place_x_cm},
                            {"tpc_place_y_cm", config.tpc_place_y_cm},
                            {"tpc_place_z_cm", config.tpc_place_z_cm},
                            {"tpc_rotation_x_rad", config.tpc_rotation_x_rad},
                            {"tpc_rotation_y_rad", config.tpc_rotation_y_rad},
                            {"tpc_rotation_z_rad", config.tpc_rotation_z_rad},
                            {"nominal_electric_field_vcm", provenance.nominal_electric_field_vcm},
                            {"space_charge_scale_side0", provenance.space_charge_scale_side0},
                            {"space_charge_scale_side1", provenance.space_charge_scale_side1},
                            {"magnetic_field_scale", provenance.magnetic_field_scale}});
  write_parameters<int>({{"zero_magnetic_field", provenance.zero_magnetic_field ? 1 : 0},
                         {"maximum_integration_steps", static_cast<int>(PHGarfieldForwardTransport::MaxForwardDriftSteps)}});
  write_parameters<double>({{"maximum_drift_time_ns", summary.maximum_drift_time_ns}});
  write_parameters<int>({{"cell_center_probes", 1}});
  write_named({{"reach_probe_trajectories", std::to_string(summary.probe_trajectories)},
               {"cdb_timestamp", std::to_string(config.cdb_timestamp)}});

  return true;
}

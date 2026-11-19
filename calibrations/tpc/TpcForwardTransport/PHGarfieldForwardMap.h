#ifndef PHGARFIELD_FORWARDMAP_H
#define PHGARFIELD_FORWARDMAP_H

#include "PHGarfieldForwardTransport.h"

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>

class TDirectory;

class PHGarfieldForwardMap
{
 public:
  static constexpr int SchemaVersion = 9;
  static constexpr const char* SchemaName =
      "PHG4TpcGarfieldAdaptiveTransportLUT";
  static constexpr const char* TransportContract =
      "garfield_full_forward_mean_transport_v1";
  static constexpr const char* SpaceChargeZCoordinatePolicy =
      "max_zero_abs_z_minus_cm_half_width_v1";
  static constexpr const char* MagneticFieldDomainPolicy =
      "cartesian_domain_full_tpc_box_and_step_checked_v1";
  static constexpr const char* NominalTimeReferencePolicy =
      "r_mid_phi0_cm_surface_each_side_v1";

  struct Config
  {
    std::string output_filename;
    std::string configuration_label;

    unsigned int phi_nodes{36};
    unsigned int radial_nodes{24};
    unsigned int longitudinal_nodes{28};

    double inner_radius_cm{21.6};
    double outer_radius_cm{76.4};
    double max_radial_field_extension_cm{0.3};  // Bounded nearest-edge extension; gas geometry is unchanged.
    double central_membrane_half_width_cm{0.28};
    double readout_abs_z_cm{102.605};
    double integration_step_ns{12.5};
    double nominal_phg4_drift_velocity_cm_per_ns{0.00755};
    double maximum_nominal_time_mismatch{0.02};  // Warning threshold for the gas-model/reference time difference.
    double interpolation_tolerance_um{0.0};  // CLI resolves zero to the adaptive guard default.
    double tpc_place_x_cm{0.0};
    double tpc_place_y_cm{0.0};
    double tpc_place_z_cm{0.0};
    double tpc_rotation_x_rad{0.0};
    double tpc_rotation_y_rad{0.0};
    double tpc_rotation_z_rad{0.0};

    std::string cdb_global_tag;
    std::uint64_t cdb_timestamp{0};
    bool overwrite{false};
  };

  struct NodeValue
  {
    bool reaches_readout{false};
    double delta_r_cm{0.0};
    double delta_rphi_cm{0.0};
    double drift_time_ns{0.0};
  };

  struct Provenance
  {
    std::string field_model;
    std::array<std::string, 2> electric_field_map_3d{};
    std::array<std::string, 2> electric_field_map_3d_md5{};
    std::string frame_field_map;
    std::string frame_field_map_md5;
    std::array<std::string, 2> frame_field_map_3d{};
    std::array<std::string, 2> frame_field_map_3d_md5{};
    std::array<double, 2> frame_charge_scale{{0.0, 0.0}};
    std::array<double, 2> cm_field_factor{{1.0, 1.0}};
    std::array<double, 4> field_cage_offsets_v{{0.0, 0.0, 0.0, 0.0}};
    bool use_ifc_voltage_distortion{false};
    bool use_ofc_voltage_distortion{false};

    double nominal_electric_field_vcm{0.0};
    double space_charge_scale_side0{0.0};
    double space_charge_scale_side1{0.0};
    double magnetic_field_scale{0.0};
    bool zero_magnetic_field{false};
    std::string electric_field_map;
    std::string magnetic_field_map;
    std::string gas_file;
    std::string electric_field_map_md5;
    std::string magnetic_field_map_md5;
    std::string gas_file_md5;
    unsigned int maximum_integration_steps{0};
    std::string reach_probe_scheme;
    double scale_zero_reference_time_side0_ns{0.0};
    double scale_zero_reference_time_side1_ns{0.0};
    double scale_zero_maximum_time_mismatch{0.0};
  };

  struct Summary
  {
    unsigned long long probe_trajectories{0};
    double maximum_drift_time_ns{0.0};
    double scale_zero_reference_time_side0_ns{0.0};
    double scale_zero_reference_time_side1_ns{0.0};
    double scale_zero_maximum_time_mismatch{0.0};
  };

  static bool EncodeNode(double r0_cm,
                         double phi0_rad,
                         const PHGarfieldForwardTransport::ForwardDriftResult& result,
                         NodeValue& value);

  static bool ValidateConfig(const Config& config, std::ostream& error);

  // Files use their MD5; gas directories hash the ordered MERGED_E###.gas
  // inventory used by PHGarfield. Empty or unreadable payloads return empty.
  static std::string PayloadMd5(const std::string& filename);
  // Stable configuration/manifest fingerprint shared by tools and readers.
  static std::string Fingerprint(const std::string& text);

  static std::string BuildCanonicalConfiguration(
      const Config& config,
      const Provenance& provenance);

  // Shared physical validation, timing reference and actual payload capture
  // used by direct adaptive LUT generation.
  static bool Prepare(const PHGarfieldForwardTransport& garfield, const Config& config,
                      Provenance& provenance, Summary& summary,
                      std::ostream& log, std::ostream& error);
  static bool WriteMetadata(TDirectory* directory, const Config& config,
                            const Provenance& provenance, const Summary& summary);
};

#endif

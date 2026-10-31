#ifndef PHGARFIELD_FORWARDMAP_H
#define PHGARFIELD_FORWARDMAP_H

#include "PHGarfield.h"

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>

class TDirectory;
class TFile;

class PHGarfieldForwardMap
{
 public:
  static constexpr int SchemaVersion = 1;           // historical axisymmetric LUTs
  static constexpr int ExtendedSchemaVersion = 2;   // explicit 3D/frame provenance
  static constexpr int NumericalSchemaVersion = 3;  // explicit interpolation-error cell mask
  static constexpr const char* SchemaName =
      "PHG4TpcGarfieldFullTransportLUT";
  static constexpr const char* TransportContract =
      "garfield_full_forward_mean_transport_v1";
  static constexpr const char* SpaceChargeZCoordinatePolicy =
      "max_zero_abs_z_minus_cm_half_width_v1";
  static constexpr const char* MagneticFieldDomainPolicy =
      "complete_cartesian_lattice_full_tpc_box_and_step_checked_v1";
  static constexpr const char* NominalTimeReferencePolicy =
      "r_mid_phi0_cm_surface_each_side_v1";

  enum class ReachCellClass
  {
    AllLost = 0,
    Mixed = 1,
    AllReach = 2,
    NumericallyUnresolved = 3
  };

  struct Config
  {
    std::string output_filename;
    std::string configuration_label;

    unsigned int phi_nodes{36};
    unsigned int radial_nodes{24};
    unsigned int longitudinal_nodes{28};
    // Execution partition only; physical axes and canonical configuration are shared.
    unsigned int phi_chunk_index{0};
    unsigned int phi_chunk_count{1};

    double inner_radius_cm{21.6};
    double outer_radius_cm{76.4};
    double central_membrane_half_width_cm{0.28};
    double readout_abs_z_cm{102.605};
    double integration_step_ns{12.5};
    double nominal_phg4_drift_velocity_cm_per_ns{0.00755};
    double maximum_nominal_time_mismatch{0.02};
    double interpolation_tolerance_um{0.0};  // zero preserves historical schemas
    double tpc_place_x_cm{0.0};
    double tpc_place_y_cm{0.0};
    double tpc_place_z_cm{0.0};
    double tpc_rotation_x_rad{0.0};
    double tpc_rotation_y_rad{0.0};
    double tpc_rotation_z_rad{0.0};

    std::string cdb_global_tag;
    std::uint64_t cdb_timestamp{0};
    bool probe_cell_centers{true};
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
    unsigned long long source_nodes{0};
    unsigned long long reaching_nodes{0};
    unsigned long long lost_nodes{0};
    unsigned long long probe_trajectories{0};
    unsigned long long all_reach_cells{0};
    unsigned long long all_lost_cells{0};
    unsigned long long mixed_cells{0};
    unsigned long long numerically_unresolved_cells{0};
    double maximum_drift_time_ns{0.0};
    double scale_zero_reference_time_side0_ns{0.0};
    double scale_zero_reference_time_side1_ns{0.0};
    double scale_zero_maximum_time_mismatch{0.0};
  };

  static bool EncodeNode(double r0_cm,
                         double phi0_rad,
                         const PHGarfield::ForwardDriftResult& result,
                         NodeValue& value);

  static bool ValidateConfig(const Config& config, std::ostream& error);
  // Shared schema1-3 axes, guards, endpoints, reachability and mask checks.
  static bool ValidateRegular(TFile& file, std::ostream& error);

  // Files use their MD5; gas directories hash the ordered MERGED_E###.gas
  // inventory used by PHGarfield. Empty or unreadable payloads return empty.
  static std::string PayloadMd5(const std::string& filename);
  // Stable configuration/manifest fingerprint shared by tools and readers.
  static std::string Fingerprint(const std::string& text);

  static std::string BuildCanonicalConfiguration(
      const Config& config,
      const Provenance& provenance,
      int schema_version = SchemaVersion);

  // Shared physical validation, timing reference and actual payload capture.
  // Schemas1-3 and adaptive generation use this same preflight definition.
  static bool Prepare(const PHGarfield&, const Config&, Provenance&, Summary&,
                      std::ostream& log, std::ostream& error);
  static bool WriteMetadata(TDirectory*, const Config&, const Provenance&,
                            const Summary&, int schema_version,
                            const std::string& schema_name = "");

  static bool Generate(const PHGarfield& garfield,
                       const Config& config,
                       Summary& summary,
                       std::ostream& log,
                       std::ostream& error);
};

#endif

#ifndef PHGARFIELDFORWARDSETUP_H
#define PHGARFIELDFORWARDSETUP_H

#include <array>
#include <functional>
#include <string>

// Internal input-selection policy for PHGarfield's opt-in forward transport.
// Resource loading, field physics, ownership, and node publication stay in
// PHGarfield. The callbacks keep this policy independent of CDB and ROOT.
namespace PHGarfieldForwardSetup
{
  using CdbLookup = std::function<std::string(const std::string&)>;
  using FileExists = std::function<bool(const std::string&)>;

  enum class MagneticFieldSource
  {
    Zero,
    ExplicitFile,
    NodeTree,
    CdbFile
  };

  struct Request
  {
    bool has_node_tree{false};
    std::array<bool, 2> space_charge_scale_override{{false, false}};
    bool zero_magnetic_field{false};
    std::string magnetic_field_override;
    bool load_tpc_fee_map{true};

    bool electric_field_map_override{false};
    bool use_2d_electric_field_map{false};
    std::array<std::string, 2> electric_field_maps_3d;
    std::array<bool, 2> electric_field_map_3d_override{{false, false}};
    bool frame_electric_field_map_override{false};
    std::array<std::string, 2> frame_electric_field_maps_3d;
    std::array<bool, 2> frame_electric_field_map_3d_override{{false, false}};
  };

  struct Selection
  {
    bool valid{false};
    bool use_run_dependent_inputs{false};
    bool load_tpc_fee_map{true};
    MagneticFieldSource magnetic_field_source{MagneticFieldSource::CdbFile};
    std::string magnetic_field_override;
    std::array<std::string, 2> electric_field_maps_3d;
    std::array<std::string, 2> frame_electric_field_maps_3d;
  };

  Selection SelectInputs(const Request& request);

  // Keep the lookups separate: magnetic-field resource loading can fail
  // before the fee lookup, and run-dependent CDB loading precedes gas lookup.
  std::string ResolveMagneticField(const Selection& selection, const CdbLookup& lookup);
  std::string ResolveFeeMap(const Selection& selection, const CdbLookup& lookup);

  struct GasRequest
  {
    std::string override_filename;
    std::string default_filename;
  };

  enum class GasSource
  {
    ExplicitFile,
    CdbFile,
    DefaultFile
  };

  struct GasSelection
  {
    bool valid{true};
    GasSource source{GasSource::CdbFile};
    std::string filename;
    // Retained for the existing missing-CDB-payload diagnostic on fallback.
    std::string unavailable_cdb_filename;
  };

  GasSelection ResolveGas(const GasRequest& request,
                          const CdbLookup& lookup,
                          const FileExists& exists);
}  // namespace PHGarfieldForwardSetup

#endif

#include "../PHGarfieldForwardSetup.h"

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
  using namespace PHGarfieldForwardSetup;

  bool Check(bool condition, const char* message)
  {
    if (!condition)
    {
      std::cerr << message << '\n';
    }
    return condition;
  }

  struct Inputs
  {
    std::map<std::string, std::string> payloads;
    std::vector<std::string> queries;
    std::vector<std::string> files_checked;
    std::vector<std::string> existing_files;

    std::string Lookup(const std::string& payload)
    {
      queries.push_back(payload);
      return payloads[payload];
    }

    bool Exists(const std::string& filename)
    {
      files_checked.push_back(filename);
      for (const auto& existing : existing_files)
      {
        if (existing == filename)
        {
          return true;
        }
      }
      return false;
    }
  };

  bool CheckStandaloneAndMagneticSelection()
  {
    const std::array<std::array<bool, 2>, 4> overrides{{{{false, false}}, {{true, false}}, {{false, true}}, {{true, true}}}};
    for (const auto& side_overrides : overrides)
    {
      Request request;
      request.space_charge_scale_override = side_overrides;
      const auto standalone = SelectInputs(request);
      const bool both = side_overrides[0] && side_overrides[1];
      if (!Check(standalone.valid == both,
                 "Standalone input selection must require explicit coefficients on both sides") ||
          !Check(!standalone.use_run_dependent_inputs,
                 "Standalone input selection must skip reconstruction CDB and run-cage configuration"))
      {
        return false;
      }

      request.has_node_tree = true;
      const auto node_selection = SelectInputs(request);
      if (!Check(node_selection.valid && node_selection.use_run_dependent_inputs,
                 "Node-tree initialization must retain run-dependent inputs even with both sides pinned"))
      {
        return false;
      }
    }

    struct MagneticCase
    {
      bool zero_field;
      std::string explicit_file;
      bool has_node_tree;
      MagneticFieldSource expected_source;
      std::string expected_path;
      bool expect_cdb;
    };
    const std::array<MagneticCase, 7> cases{{{true, "pinned.root", true, MagneticFieldSource::Zero, "ZERO_UNIFORM_FIELD", false},
                                             {true, "", false, MagneticFieldSource::Zero, "ZERO_UNIFORM_FIELD", false},
                                             {false, "pinned.root", true, MagneticFieldSource::ExplicitFile, "pinned.root", false},
                                             {false, "pinned.root", false, MagneticFieldSource::ExplicitFile, "pinned.root", false},
                                             {false, "", true, MagneticFieldSource::NodeTree, "FIELD_MAP node", false},
                                             {false, "", false, MagneticFieldSource::CdbFile, "tracking.root", true},
                                             {true, "pinned.root", false, MagneticFieldSource::Zero, "ZERO_UNIFORM_FIELD", false}}};
    for (const auto& test : cases)
    {
      Request request;
      request.has_node_tree = test.has_node_tree;
      request.space_charge_scale_override = {{true, true}};
      request.zero_magnetic_field = test.zero_field;
      request.magnetic_field_override = test.explicit_file;
      Inputs inputs;
      inputs.payloads["FIELDMAP_TRACKING"] = "tracking.root";
      const auto selection = SelectInputs(request);
      const auto resolved = ResolveMagneticField(selection, [&inputs](const std::string& payload)
                                                 { return inputs.Lookup(payload); });
      const std::vector<std::string> expected_queries = test.expect_cdb
                                                            ? std::vector<std::string>{"FIELDMAP_TRACKING"}
                                                            : std::vector<std::string>{};
      if (!Check(selection.magnetic_field_source == test.expected_source &&
                     resolved == test.expected_path && inputs.queries == expected_queries,
                 "Magnetic-field priority or CDB isolation changed"))
      {
        return false;
      }
    }
    return true;
  }

  bool CheckElectricFieldFamilies()
  {
    Request request;
    request.has_node_tree = true;
    request.electric_field_maps_3d = {{"space0.root", "space1.root"}};
    request.frame_electric_field_maps_3d = {{"frame0.root", "frame1.root"}};
    auto selection = SelectInputs(request);
    if (!Check(selection.electric_field_maps_3d == request.electric_field_maps_3d &&
                   selection.frame_electric_field_maps_3d == request.frame_electric_field_maps_3d,
               "Automatic 3D selections must survive when no 2D family is selected"))
    {
      return false;
    }

    request.electric_field_map_override = true;
    selection = SelectInputs(request);
    const std::array<std::string, 2> empty{{"", ""}};
    if (!Check(selection.electric_field_maps_3d == empty &&
                   selection.frame_electric_field_maps_3d == empty,
               "An explicit 2D electric map must discard stale automatic 3D electric and frame maps"))
    {
      return false;
    }

    request.electric_field_map_3d_override[0] = true;
    request.frame_electric_field_map_3d_override[1] = true;
    selection = SelectInputs(request);
    if (!Check(selection.electric_field_maps_3d == std::array<std::string, 2>{{"space0.root", ""}} &&
                   selection.frame_electric_field_maps_3d == request.frame_electric_field_maps_3d,
               "Pinned side maps must retain per-side precedence and the other automatic frame side"))
    {
      return false;
    }

    request.frame_electric_field_map_override = true;
    selection = SelectInputs(request);
    if (!Check(selection.frame_electric_field_maps_3d == std::array<std::string, 2>{{"", "frame1.root"}},
               "An explicit 2D frame map must discard only unpinned 3D frame paths"))
    {
      return false;
    }

    request.electric_field_map_override = false;
    request.frame_electric_field_map_override = false;
    request.use_2d_electric_field_map = true;
    selection = SelectInputs(request);
    if (!Check(selection.electric_field_maps_3d == std::array<std::string, 2>{{"space0.root", ""}} &&
                   selection.frame_electric_field_maps_3d == request.frame_electric_field_maps_3d,
               "CDB 2D electric selection must preserve independent frame-map selection"))
    {
      return false;
    }

    request.electric_field_maps_3d = selection.electric_field_maps_3d;
    request.frame_electric_field_maps_3d = selection.frame_electric_field_maps_3d;
    const auto repeated = SelectInputs(request);
    return Check(repeated.electric_field_maps_3d == selection.electric_field_maps_3d &&
                     repeated.frame_electric_field_maps_3d == selection.frame_electric_field_maps_3d,
                 "Repeated selection must retain pinned paths and keep discarded paths empty");
  }

  bool CheckFeeAndGasResolution()
  {
    Request request;
    request.has_node_tree = true;
    Inputs inputs;
    inputs.payloads["TPC_FEE_CHANNEL_MAP"] = "fee.root";
    const auto lookup = [&inputs](const std::string& payload)
    { return inputs.Lookup(payload); };
    const auto exists = [&inputs](const std::string& filename)
    { return inputs.Exists(filename); };

    request.load_tpc_fee_map = false;
    if (!Check(ResolveFeeMap(SelectInputs(request), lookup).empty() && inputs.queries.empty(),
               "Disabling the plotting fee map must eliminate its CDB lookup"))
    {
      return false;
    }
    request.load_tpc_fee_map = true;
    if (!Check(ResolveFeeMap(SelectInputs(request), lookup) == "fee.root" &&
                   inputs.queries == std::vector<std::string>{"TPC_FEE_CHANNEL_MAP"},
               "Enabled fee-map loading must retain its existing payload key"))
    {
      return false;
    }

    inputs.queries.clear();
    GasRequest gas{"pinned.gas", "default.gas"};
    inputs.existing_files = {"pinned.gas"};
    auto resolved = ResolveGas(gas, lookup, exists);
    if (!Check(resolved.valid && resolved.source == GasSource::ExplicitFile &&
                   resolved.filename == "pinned.gas" && inputs.queries.empty(),
               "An existing explicit gas file must avoid gas CDB lookup"))
    {
      return false;
    }

    inputs.files_checked.clear();
    inputs.existing_files.clear();
    resolved = ResolveGas(gas, lookup, exists);
    if (!Check(!resolved.valid && resolved.filename == "pinned.gas" &&
                   resolved.source == GasSource::ExplicitFile && inputs.queries.empty() &&
                   inputs.files_checked == std::vector<std::string>{"pinned.gas"},
               "A missing explicit gas file must fail without CDB or default substitution"))
    {
      return false;
    }

    gas.override_filename.clear();
    inputs.payloads["PHGARFIELD_GAS"] = "cdb.gas";
    inputs.existing_files = {"cdb.gas"};
    inputs.files_checked.clear();
    resolved = ResolveGas(gas, lookup, exists);
    if (!Check(resolved.valid && resolved.source == GasSource::CdbFile &&
                   resolved.filename == "cdb.gas" &&
                   inputs.queries == std::vector<std::string>{"PHGARFIELD_GAS"},
               "An available CDB gas file must be selected with the existing payload key"))
    {
      return false;
    }

    inputs.existing_files.clear();
    inputs.files_checked.clear();
    resolved = ResolveGas(gas, lookup, exists);
    if (!Check(resolved.valid && resolved.source == GasSource::DefaultFile &&
                   resolved.filename == "default.gas" && resolved.unavailable_cdb_filename == "cdb.gas" &&
                   inputs.files_checked == std::vector<std::string>{"cdb.gas"},
               "An unavailable CDB gas file must fall back without validating the default file"))
    {
      return false;
    }

    inputs.payloads["PHGARFIELD_GAS"].clear();
    inputs.files_checked.clear();
    resolved = ResolveGas(gas, lookup, exists);
    return Check(resolved.valid && resolved.source == GasSource::DefaultFile &&
                     resolved.filename == "default.gas" && resolved.unavailable_cdb_filename.empty() &&
                     inputs.files_checked.empty(),
                 "An absent CDB gas payload must use the default and leave loading errors to PHGarfield");
  }
}  // namespace

int main()
{
  return CheckStandaloneAndMagneticSelection() &&
                 CheckElectricFieldFamilies() &&
                 CheckFeeAndGasResolution()
             ? 0
             : 1;
}

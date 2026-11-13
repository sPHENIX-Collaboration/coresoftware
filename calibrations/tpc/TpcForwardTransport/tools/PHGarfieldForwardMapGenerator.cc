#include "PHGarfieldAdaptiveMap.h"
#include "PHGarfieldForwardMap.h"
#include "PHGarfieldForwardTransport.h"
#include "PHGarfieldMapTools.h"

#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace
{
  struct Options
  {
    PHGarfieldForwardMap::Config map;
    PHGarfieldAdaptiveMap::Policy adaptive_policy;
    std::uint64_t adaptive_first{0}, adaptive_end{0};

    std::array<std::string, 2> electric_field_map_3d{};
    std::array<std::string, 2> frame_field_map_3d{};
    double frame_charge_scale{0.0};
    std::string gas_file;
    std::string magnetic_field_map;
    std::string global_tag{"FieldMapTest"};
    std::uint64_t timestamp{1};
    double magnetic_field_scale{1.0};
    double space_charge_scale_side0{1.0};
    double space_charge_scale_side1{1.0};
    double nominal_electric_field_vcm{400.0};
    bool zero_magnetic_field{false};
    bool help_requested{false};
    bool check_inputs{false};
    int verbosity{0};
  };

  using Destination = std::variant<bool*, std::string*, double*, unsigned int*, std::uint64_t*, int*>;
  struct Option
  {
    const char* name;
    Destination destination;
    const char* argument;
    const char* help;
    const char* section{};
  };

  // Keep accepted options, destination types, and help in one place.
  std::vector<Option> option_table(Options& o)
  {
    return {
        {"--electric-field-map-3d-side0", o.electric_field_map_3d.data(), "FILE", "Override the south 3D space-charge CDB map", "Physics inputs"},
        {"--electric-field-map-3d-side1", &o.electric_field_map_3d[1], "FILE", "Override the north 3D space-charge CDB map"},
        {"--frame-field-map-3d-side0", o.frame_field_map_3d.data(), "FILE", "Override the south 3D frame CDB map"},
        {"--frame-field-map-3d-side1", &o.frame_field_map_3d[1], "FILE", "Override the north 3D frame CDB map"},
        {"--frame-charge-scale", &o.frame_charge_scale, "VALUE", "Frame coefficient (default 0)"},
        {"--gas-file", &o.gas_file, "FILE", "Override the PHGARFIELD_GAS CDB payload"},
        {"--magnetic-field-map", &o.magnetic_field_map, "FILE", "Override the FIELDMAP_TRACKING CDB payload"},
        {"--magnetic-field-scale", &o.magnetic_field_scale, "VALUE", "Magnetic-field map scale (default 1)"},
        {"--zero-b-field", &o.zero_magnetic_field, "", "Generate with B=0"},
        {"--scale-side0", &o.space_charge_scale_side0, "VALUE", "Negative-z correction scale (default 1)"},
        {"--scale-side1", &o.space_charge_scale_side1, "VALUE", "Positive-z correction scale (default 1)"},
        {"--nominal-e-vcm", &o.nominal_electric_field_vcm, "VALUE", "Nominal longitudinal E magnitude (default 400)"},
        {"--step-ns", &o.map.integration_step_ns, "VALUE", "Euler step in ns (default 12.5)", "Grid and integration"},
        {"--phi-nodes", &o.map.phi_nodes, "N", "Physical phi nodes (default 36)"},
        {"--radial-nodes", &o.map.radial_nodes, "N", "Physical radial nodes (default 24)"},
        {"--longitudinal-nodes", &o.map.longitudinal_nodes, "N", "Physical z nodes per side (default 28)"},
        {"--inner-radius-cm", &o.map.inner_radius_cm, "VALUE", "Gas inner radius (default 21.6)"},
        {"--outer-radius-cm", &o.map.outer_radius_cm, "VALUE", "Gas outer radius (default 76.4)"},
        {"--radial-field-extension-cm", &o.map.max_radial_field_extension_cm, "VALUE", "Maximum nearest-edge radial field extension (default 0.3; 0 is strict)"},
        {"--cm-half-width-cm", &o.map.central_membrane_half_width_cm, "VALUE", "Central-membrane half-width (default 0.28)"},
        {"--readout-abs-z-cm", &o.map.readout_abs_z_cm, "VALUE", "Readout plane |z| (default 102.605)"},
        {"--phg4-drift-velocity", &o.map.nominal_phg4_drift_velocity_cm_per_ns, "VALUE", "PHG4 calibration velocity in cm/ns (default 0.00755)"},
        {"--max-nominal-time-mismatch", &o.map.maximum_nominal_time_mismatch, "F", "Reference drift-time warning threshold (fraction; default 0.02)"},
        {"--tpc-place-x-cm", &o.map.tpc_place_x_cm, "VALUE", "TPC global x placement (default 0)"},
        {"--tpc-place-y-cm", &o.map.tpc_place_y_cm, "VALUE", "TPC global y placement (default 0)"},
        {"--tpc-place-z-cm", &o.map.tpc_place_z_cm, "VALUE", "TPC global z placement (default 0)"},
        {"--tpc-rotation-x-rad", &o.map.tpc_rotation_x_rad, "VALUE", "TPC x rotation (default 0)"},
        {"--tpc-rotation-y-rad", &o.map.tpc_rotation_y_rad, "VALUE", "TPC y rotation (default 0)"},
        {"--tpc-rotation-z-rad", &o.map.tpc_rotation_z_rad, "VALUE", "TPC z rotation (default 0)"},
        {"--adaptive-root-first", &o.adaptive_first, "N", "Inclusive coarse root ID (default 0)", "Advanced refinement and resource controls"},
        {"--adaptive-root-end", &o.adaptive_end, "N", "Exclusive root ID (default all roots)"},
        {"--adaptive-max-depth", &o.adaptive_policy.max_depth, "N", "Per-root depth cap (default 12)"},
        {"--adaptive-max-nodes-per-root", &o.adaptive_policy.max_nodes_per_root, "N", "Per-root node cap (default 127)"},
        {"--adaptive-max-queries-per-root", &o.adaptive_policy.max_queries_per_root, "N", "Unique stencil-query cap (default 2048)"},
        {"--adaptive-max-total-points", &o.adaptive_policy.max_total_points, "N", "Global per-file point cap (default 5000000)"},
        {"--adaptive-max-total-cells", &o.adaptive_policy.max_total_cells, "N", "Global per-file cell cap (default 2000000)"},
        {"--adaptive-max-bytes", &o.adaptive_policy.max_estimated_bytes, "N", "Estimated construction-allocation cap (default 2147483648)"},
        {"--interpolation-tolerance-um", &o.map.interpolation_tolerance_um, "V", "Maximum actual 27-stencil error (adaptive default 10)"},
        {"--output", &o.map.output_filename, "FILE", "Output ROOT file (also -o)", "CDB and output"},
        {"--global-tag", &o.global_tag, "TAG", "CDB global tag (default FieldMapTest)"},
        {"--timestamp", &o.timestamp, "VALUE", "CDB timestamp (default 1)"},
        {"--label", &o.map.configuration_label, "TEXT", "Human-readable configuration label"},
        {"--force", &o.map.overwrite, "", "Overwrite an existing output file"},
        {"--verbosity", &o.verbosity, "N", "PHGarfield verbosity"},
        {"--check-inputs", &o.check_inputs, "", "Validate inputs and timing without generating a LUT"},
        {"--help", &o.help_requested, "", "Show this message (also -h)"}};
  }

  void print_usage(std::ostream& output, const char* program)
  {
    output << "Usage: " << program << " --output FILE [--global-tag TAG --timestamp VALUE] [options]\n"
           << "\nGenerate a ready-to-use adaptive PHGarfield LUT directly; no conversion step.\n"
           << "Missing map paths use CDB. Frame maps are requested by a nonzero frame coefficient or an explicit frame path.\n";
    Options defaults;
    const auto flags = output.flags();
    for (const auto& option : option_table(defaults))
    {
      if (!option.help[0])
      {
        continue;
      }
      if (option.section)
      {
        output << '\n'
               << option.section << ":\n";
      }
      std::string label = option.name;
      if (option.argument[0])
      {
        label += std::string(" ") + option.argument;
      }
      output << "  " << std::left << std::setw(40) << label << option.help << '\n';
    }
    output.flags(flags);
  }

  void assign_option(const Option& option, const std::string& text)
  {
    std::visit([&](auto destination)
               {
      using Pointer = decltype(destination);
      using Value = std::remove_pointer_t<Pointer>;
      if constexpr (std::is_same_v<Value, bool>)
      {
        *destination = true;
      }
      else if constexpr (std::is_same_v<Value, std::string>)
      {
        *destination = text;
      }
      else
      {
        std::size_t parsed = 0;
        const auto invalid = [&](const char* type)
        {
          throw std::runtime_error(std::string("Invalid ") + type + " value after " + option.name);
        };
        if constexpr (std::is_unsigned_v<Value>)
        {
          if (text.empty() || text.front() == '-')
          {
            invalid("unsigned");
          }
          const auto value = std::stoull(text, &parsed);
          if (parsed != text.size() || value > std::numeric_limits<Value>::max())
          {
            invalid("unsigned");
          }
          *destination = static_cast<Value>(value);
        }
        else if constexpr (std::is_same_v<Value, int>)
        {
          const auto value = std::stoll(text, &parsed);
          if (parsed != text.size() || value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
          {
            invalid("integer");
          }
          *destination = static_cast<int>(value);
        }
        else
        {
          const auto value = std::stod(text, &parsed);
          if (parsed != text.size())
          {
            invalid("floating-point");
          }
          *destination = value;
        }
      } }, option.destination);
  }

  bool parse_options(const int argc, char* argv[], Options& options, std::ostream& error)
  {
    const auto table = option_table(options);
    try
    {
      for (int index = 1; index < argc; ++index)
      {
        const std::string argument = argv[index];
        std::string name = argument;
        if (argument == "-o")
        {
          name = "--output";
        }
        else if (argument == "-h")
        {
          name = "--help";
        }
        const auto option = std::find_if(table.begin(), table.end(), [&](const Option& entry)
                                         { return name == entry.name; });
        if (option == table.end())
        {
          throw std::runtime_error("Unknown option " + argument);
        }
        std::string text;
        if (option->argument[0])
        {
          if (index + 1 >= argc)
          {
            throw std::runtime_error("Missing value after " + argument);
          }
          text = argv[++index];
        }
        assign_option(*option, text);
        if (options.help_requested)
        {
          print_usage(std::cout, argv[0]);
          return true;
        }
      }
    }
    catch (const std::exception& exception)
    {
      error << exception.what() << "\n";
      return false;
    }

    if (options.map.interpolation_tolerance_um == 0.)
    {
      options.map.interpolation_tolerance_um = options.adaptive_policy.guard_um;
    }
    options.adaptive_policy.guard_um = options.map.interpolation_tolerance_um;
    const auto roots = PHGarfieldAdaptiveMap::RootCount(PHGarfieldAdaptiveMap::FromConfig(options.map));
    if (!options.adaptive_end)
    {
      options.adaptive_end = roots;
    }
    if (options.adaptive_first >= options.adaptive_end || options.adaptive_end > roots)
    {
      error << "Adaptive generation requires a nonempty valid root range\n";
      return false;
    }
    if (options.map.output_filename.empty())
    {
      error << "--output is required\n";
      return false;
    }
    return true;
  }

}  // namespace

int main(const int argc, char* argv[])
try
{
  Options options;
  if (!parse_options(argc, argv, options, std::cerr))
  {
    if (argc <= 1)
    {
      print_usage(std::cerr, argv[0]);
    }
    return 2;
  }
  if (options.help_requested)
  {
    return 0;
  }

  const auto require_regular_file =
      [](const std::string& filename, const char* option)
  {
    if (filename.empty())
    {
      return true;
    }
    std::error_code error;
    if (!std::filesystem::is_regular_file(filename, error) || error)
    {
      std::cerr << option << " is not a readable regular file: "
                << filename << "\n";
      return false;
    }
    return true;
  };
  if (!require_regular_file(options.magnetic_field_map, "--magnetic-field-map"))
  {
    return 2;
  }
  for (const auto& paths : {options.electric_field_map_3d, options.frame_field_map_3d})
  {
    for (const auto& path : paths)
    {
      if (!require_regular_file(path, "3D field payload"))
      {
        return 2;
      }
    }
  }
  if (!options.gas_file.empty())
  {
    std::error_code error;
    const bool valid_gas_path =
        std::filesystem::is_regular_file(options.gas_file, error) ||
        std::filesystem::is_directory(options.gas_file, error);
    if (error || !valid_gas_path)
    {
      std::cerr << "--gas-file is not a readable file or directory: "
                << options.gas_file << "\n";
      return 2;
    }
  }
  if (!std::isfinite(options.frame_charge_scale) ||
      !std::isfinite(options.magnetic_field_scale) ||
      !std::isfinite(options.space_charge_scale_side0) ||
      !std::isfinite(options.space_charge_scale_side1) ||
      !std::isfinite(options.nominal_electric_field_vcm) ||
      options.nominal_electric_field_vcm <= 0.0)
  {
    std::cerr
        << "Magnetic/space-charge scales must be finite and nominal E "
           "must be finite and positive\n";
    return 2;
  }
  if (!options.zero_magnetic_field &&
      options.magnetic_field_scale == 0.0)
  {
    std::cerr
        << "Use --zero-b-field, rather than a zero magnetic-field scale, "
           "so the map's transport contract records B=0 explicitly\n";
    return 2;
  }

  recoConsts* constants = recoConsts::instance();
  constants->set_StringFlag("CDB_GLOBALTAG", options.global_tag);
  constants->set_uint64Flag("TIMESTAMP", options.timestamp);
  options.map.cdb_global_tag = options.global_tag;
  options.map.cdb_timestamp = options.timestamp;

  std::cout << "CDB tag=" << options.global_tag << " timestamp=" << options.timestamp << '\n';
  const auto lookup = [](const std::string& payload)
  {
    const auto path = CDBInterface::instance()->getUrl(payload);
    std::cout << "CDB " << payload << ": " << path << '\n';
    return path;
  };
  options.electric_field_map_3d = PHGarfieldMapTools::resolve_field_maps(
      options.electric_field_map_3d, PHGarfieldMapTools::ElectricFieldMapPayloads, lookup);
  if (options.frame_charge_scale != 0. || !options.frame_field_map_3d[0].empty() || !options.frame_field_map_3d[1].empty())
  {
    options.frame_field_map_3d = PHGarfieldMapTools::resolve_field_maps(
        options.frame_field_map_3d, PHGarfieldMapTools::FrameFieldMapPayloads, lookup);
  }

  auto garfield = std::make_unique<PHGarfieldForwardTransport>(
      "PHGarfieldForwardMapGenerator",
      "",
      options.space_charge_scale_side0,
      options.space_charge_scale_side1,
      options.electric_field_map_3d[0], options.electric_field_map_3d[1]);
  garfield->SetForwardTransportMode(true);
  garfield->Verbosity(options.verbosity);
  garfield->SetElectricFieldMap3D(options.electric_field_map_3d[0], options.electric_field_map_3d[1]);
  garfield->SetFrameElectricFieldMap("");
  garfield->SetFrameElectricFieldMap3D(options.frame_field_map_3d[0], options.frame_field_map_3d[1]);
  garfield->SetFrameChargeScale(options.frame_charge_scale);
  garfield->SetUseSurveyGeometry(false);
  garfield->SetUseBCOkEffs(false);
  garfield->SetUsePTCorrection(false);
  garfield->SetUseCDBDriftCorrections(false);
  garfield->SetCMFieldRunFactor(1.0, 1.0);
  garfield->SetDeltaT0Tbins(0.0, 0.0);
  garfield->SetUseIFCVoltageDistortion(false);
  garfield->SetUseOFCVoltageDistortion(false);
  garfield->SetFieldCageVoltageOffsets(0.0, 0.0, 0.0, 0.0);
  garfield->SetCMVoltageDefault(options.nominal_electric_field_vcm);
  garfield->SetLoadTpcFeeMap(false);
  garfield->SetCentralMembraneHalfWidth(
      options.map.central_membrane_half_width_cm);
  garfield->SetDriftVolumeRadii(
      options.map.inner_radius_cm,
      options.map.outer_radius_cm);
  garfield->SetZeroMagneticField(options.zero_magnetic_field);
  garfield->MoveTpc(
      options.map.tpc_place_x_cm,
      options.map.tpc_place_y_cm,
      options.map.tpc_place_z_cm);
  garfield->SetTpcRotation(
      options.map.tpc_rotation_x_rad,
      options.map.tpc_rotation_y_rad,
      options.map.tpc_rotation_z_rad);
  if (!options.gas_file.empty())
  {
    garfield->SetGasFile(options.gas_file);
  }
  garfield->SetMagneticFieldMap(
      options.magnetic_field_map,
      options.magnetic_field_scale);

  if (garfield->InitRun(nullptr) != 0 || !garfield->IsReady())
  {
    std::cerr << "PHGarfield initialization failed\n";
    return 1;
  }
  for (std::size_t side = 0; side < 2; ++side)
  {
    if ((!options.electric_field_map_3d[side].empty() && !garfield->HasLoadedElectricFieldMap3D(side)) ||
        (!options.frame_field_map_3d[side].empty() && !garfield->HasLoadedFrameElectricFieldMap3D(side)))
    {
      std::cerr << "A requested 3D field payload was not loaded\n";
      return 1;
    }
  }
  PHGarfieldAdaptiveMap::BuildSummary summary;
  if (options.check_inputs)
  {
    PHGarfieldForwardMap::Provenance provenance;
    PHGarfieldForwardMap::Summary preparation;
    return PHGarfieldForwardMap::Prepare(*garfield, options.map, provenance, preparation, std::cout, std::cerr) ? 0 : 1;
  }
  return PHGarfieldAdaptiveMap::Generate(*garfield, options.map, options.adaptive_policy,
                                         options.adaptive_first, options.adaptive_end, summary, std::cout, std::cerr)
             ? 0
             : 1;
}
catch (const std::exception& error)
{
  std::cerr << "Forward-map generation failed: " << error.what() << '\n';
  return 1;
}

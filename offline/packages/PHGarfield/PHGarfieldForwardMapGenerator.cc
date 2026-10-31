#include "PHGarfield.h"
#include "PHGarfieldAdaptiveMap.h"
#include "PHGarfieldForwardMap.h"

#include <phool/recoConsts.h>

#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace
{
  struct Options
  {
    PHGarfieldForwardMap::Config map;
    PHGarfieldAdaptiveMap::Policy adaptive_policy;
    bool adaptive{false};
    std::uint64_t adaptive_first{0}, adaptive_end{0};
    bool adaptive_range_set{false};

    std::string electric_field_map;
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
    bool allow_no_electric_corrections{false};
    bool help_requested{false};
    int verbosity{0};
  };

  void print_usage(std::ostream& output, const char* program)
  {
    output
        << "Usage: " << program << " --output FILE --electric-field-map FILE [options]\n"
        << "\nGenerate a full deterministic PHGarfield forward-transport LUT.\n"
        << "\nPhysics inputs:\n"
        << "  --electric-field-map FILE       ROOT Er/Ez correction map\n"
        << "  --electric-field-map-3d-side0 FILE  Pinned south 3D space-charge map\n"
        << "  --electric-field-map-3d-side1 FILE  Pinned north 3D space-charge map\n"
        << "  --frame-field-map-3d-side0 FILE     Pinned south 3D frame map\n"
        << "  --frame-field-map-3d-side1 FILE     Pinned north 3D frame map\n"
        << "  --frame-charge-scale VALUE         Frame coefficient (default 0)\n"
        << "  --gas-file FILE                 Override the PHGARFIELD_GAS CDB payload\n"
        << "  --magnetic-field-map FILE       Override the FIELDMAP_TRACKING CDB payload\n"
        << "  --magnetic-field-scale VALUE    Magnetic-field map scale (default 1)\n"
        << "  --zero-b-field                  Generate with B=0\n"
        << "  --scale-side0 VALUE             Negative-z correction scale (default 1)\n"
        << "  --scale-side1 VALUE             Positive-z correction scale (default 1)\n"
        << "  --nominal-e-vcm VALUE           Nominal longitudinal E magnitude (default 400)\n"
        << "  --allow-no-electric-corrections Permit a nominal E/B-only map (requires both scales to be zero)\n"
        << "\nGrid and integration:\n"
        << "  --step-ns VALUE                 Euler step in ns (default 12.5)\n"
        << "  --phi-nodes N                   Physical phi nodes (default 36)\n"
        << "  --phi-chunk-index I            Execution shard index (default 0)\n"
        << "  --phi-chunk-count N            Shards; partial maps are NOT runtime-readable\n"
        << "  --radial-nodes N                Physical radial nodes (default 24)\n"
        << "  --longitudinal-nodes N          Physical z nodes per side (default 28)\n"
        << "  --inner-radius-cm VALUE         Gas inner radius (default 21.6)\n"
        << "  --outer-radius-cm VALUE         Gas outer radius (default 76.4)\n"
        << "  --cm-half-width-cm VALUE        Central-membrane half-width (default 0.28)\n"
        << "  --readout-abs-z-cm VALUE        Readout plane |z| (default 102.605)\n"
        << "  --phg4-drift-velocity VALUE     PHG4 calibration velocity in cm/ns (default 0.00755)\n"
        << "  --max-nominal-time-mismatch F   Allowed scale-zero Garfield/PHG4 fractional timing mismatch (default 0.02)\n"
        << "  --tpc-place-x-cm VALUE          TPC global x placement (default 0)\n"
        << "  --tpc-place-y-cm VALUE          TPC global y placement (default 0)\n"
        << "  --tpc-place-z-cm VALUE          TPC global z placement (default 0)\n"
        << "  --tpc-rotation-x-rad VALUE      TPC x rotation (default 0)\n"
        << "  --tpc-rotation-y-rad VALUE      TPC y rotation (default 0)\n"
        << "  --tpc-rotation-z-rad VALUE      TPC z rotation (default 0)\n"
        << "\nAdaptive anisotropic schema4 (opt-in):\n"
        << "  --adaptive                     Build a dyadic binary cell tree\n"
        << "  --adaptive-root-first N        Inclusive coarse root ID (default0)\n"
        << "  --adaptive-root-end N          Exclusive root ID (default all roots)\n"
        << "  --adaptive-max-depth N         Per-root depth cap (default12)\n"
        << "  --adaptive-max-nodes-per-root N Per-root node cap (default127)\n"
        << "  --adaptive-max-queries-per-root N Unique stencil-query cap (default2048)\n"
        << "  --adaptive-max-total-points N  Global per-file point cap (default5000000)\n"
        << "  --adaptive-max-total-cells N   Global per-file cell cap (default2000000)\n"
        << "  --adaptive-max-bytes N         Conservative global working-memory cap (default2147483648)\n"
        << "  --interpolation-tolerance-um V Maximum actual 27-stencil error (adaptive default10)\n"
        << "\nCDB and output:\n"
        << "  --global-tag TAG                CDB global tag (default FieldMapTest)\n"
        << "  --timestamp VALUE               CDB timestamp (default 1)\n"
        << "  --label TEXT                    Human-readable configuration label\n"
        << "  --force                         Overwrite an existing output file\n"
        << "  --verbosity N                   PHGarfield verbosity\n"
        << "  --help                          Show this message\n";
  }

  bool parse_options(
      const int argc,
      char* argv[],
      Options& options,
      std::ostream& error)
  {
    const auto value_after =
        [&](int& index, const std::string& option) -> std::string
    {
      if (index + 1 >= argc)
      {
        throw std::runtime_error("Missing value after " + option);
      }
      return argv[++index];
    };
    const auto unsigned_after =
        [&](int& index, const std::string& option) -> unsigned int
    {
      const std::string text = value_after(index, option);
      if (text.empty() || text.front() == '-')
      {
        throw std::runtime_error("Invalid unsigned value after " + option);
      }
      std::size_t parsed = 0;
      const unsigned long long value = std::stoull(text, &parsed);
      if (parsed != text.size() ||
          value > std::numeric_limits<unsigned int>::max())
      {
        throw std::runtime_error("Invalid unsigned value after " + option);
      }
      return static_cast<unsigned int>(value);
    };
    const auto double_after =
        [&](int& index, const std::string& option)
    {
      const std::string text = value_after(index, option);
      std::size_t parsed = 0;
      const double value = std::stod(text, &parsed);
      if (parsed != text.size())
      {
        throw std::runtime_error("Invalid floating-point value after " + option);
      }
      return value;
    };
    const auto uint64_after =
        [&](int& index, const std::string& option)
    {
      const std::string text = value_after(index, option);
      if (text.empty() || text.front() == '-')
      {
        throw std::runtime_error("Invalid unsigned value after " + option);
      }
      std::size_t parsed = 0;
      const std::uint64_t value = std::stoull(text, &parsed);
      if (parsed != text.size())
      {
        throw std::runtime_error("Invalid unsigned value after " + option);
      }
      return value;
    };
    const auto int_after =
        [&](int& index, const std::string& option)
    {
      const std::string text = value_after(index, option);
      std::size_t parsed = 0;
      const long long value = std::stoll(text, &parsed);
      if (parsed != text.size() ||
          value < std::numeric_limits<int>::min() ||
          value > std::numeric_limits<int>::max())
      {
        throw std::runtime_error("Invalid integer value after " + option);
      }
      return static_cast<int>(value);
    };

    try
    {
      for (int index = 1; index < argc; ++index)
      {
        const std::string option = argv[index];
        if (option == "--help" || option == "-h")
        {
          print_usage(std::cout, argv[0]);
          options.help_requested = true;
          return true;
        }
        if (option == "--output" || option == "-o")
        {
          options.map.output_filename = value_after(index, option);
        }
        else if (option == "--electric-field-map")
        {
          options.electric_field_map = value_after(index, option);
        }
        else if (option == "--electric-field-map-3d-side0" || option == "--electric-field-map-3d-side1")
        {
          options.electric_field_map_3d[option.back() == '0' ? 0 : 1] = value_after(index, option);
        }
        else if (option == "--frame-field-map-3d-side0" || option == "--frame-field-map-3d-side1")
        {
          options.frame_field_map_3d[option.back() == '0' ? 0 : 1] = value_after(index, option);
        }
        else if (option == "--frame-charge-scale")
        {
          options.frame_charge_scale = double_after(index, option);
        }
        else if (option == "--gas-file")
        {
          options.gas_file = value_after(index, option);
        }
        else if (option == "--magnetic-field-map")
        {
          options.magnetic_field_map = value_after(index, option);
        }
        else if (option == "--magnetic-field-scale")
        {
          options.magnetic_field_scale =
              double_after(index, option);
        }
        else if (option == "--zero-b-field")
        {
          options.zero_magnetic_field = true;
        }
        else if (option == "--scale-side0")
        {
          options.space_charge_scale_side0 =
              double_after(index, option);
        }
        else if (option == "--scale-side1")
        {
          options.space_charge_scale_side1 =
              double_after(index, option);
        }
        else if (option == "--nominal-e-vcm")
        {
          options.nominal_electric_field_vcm =
              double_after(index, option);
        }
        else if (option == "--allow-no-electric-corrections")
        {
          options.allow_no_electric_corrections = true;
        }
        else if (option == "--step-ns")
        {
          options.map.integration_step_ns =
              double_after(index, option);
        }
        else if (option == "--adaptive")
        {
          options.adaptive = true;
        }
        else if (option == "--adaptive-root-first")
        {
          options.adaptive_first = uint64_after(index, option);
          options.adaptive_range_set = true;
        }
        else if (option == "--adaptive-root-end")
        {
          options.adaptive_end = uint64_after(index, option);
          options.adaptive_range_set = true;
        }
        else if (option == "--adaptive-max-depth")
        {
          options.adaptive_policy.max_depth = unsigned_after(index, option);
        }
        else if (option == "--adaptive-max-nodes-per-root")
        {
          options.adaptive_policy.max_nodes_per_root = unsigned_after(index, option);
        }
        else if (option == "--adaptive-max-queries-per-root")
        {
          options.adaptive_policy.max_queries_per_root = unsigned_after(index, option);
        }
        else if (option == "--adaptive-max-total-points")
        {
          options.adaptive_policy.max_total_points = uint64_after(index, option);
        }
        else if (option == "--adaptive-max-total-cells")
        {
          options.adaptive_policy.max_total_cells = uint64_after(index, option);
        }
        else if (option == "--adaptive-max-bytes")
        {
          options.adaptive_policy.max_estimated_bytes = uint64_after(index, option);
        }
        else if (option == "--interpolation-tolerance-um")
        {
          options.map.interpolation_tolerance_um = double_after(index, option);
        }
        else if (option == "--phi-chunk-index")
        {
          options.map.phi_chunk_index = unsigned_after(index, option);
        }
        else if (option == "--phi-chunk-count")
        {
          options.map.phi_chunk_count = unsigned_after(index, option);
        }
        else if (option == "--phi-nodes")
        {
          options.map.phi_nodes = unsigned_after(index, option);
        }
        else if (option == "--radial-nodes")
        {
          options.map.radial_nodes = unsigned_after(index, option);
        }
        else if (option == "--longitudinal-nodes")
        {
          options.map.longitudinal_nodes = unsigned_after(index, option);
        }
        else if (option == "--inner-radius-cm")
        {
          options.map.inner_radius_cm =
              double_after(index, option);
        }
        else if (option == "--outer-radius-cm")
        {
          options.map.outer_radius_cm =
              double_after(index, option);
        }
        else if (option == "--cm-half-width-cm")
        {
          options.map.central_membrane_half_width_cm =
              double_after(index, option);
        }
        else if (option == "--readout-abs-z-cm")
        {
          options.map.readout_abs_z_cm =
              double_after(index, option);
        }
        else if (option == "--phg4-drift-velocity")
        {
          options.map.nominal_phg4_drift_velocity_cm_per_ns =
              double_after(index, option);
        }
        else if (option == "--max-nominal-time-mismatch")
        {
          options.map.maximum_nominal_time_mismatch =
              double_after(index, option);
        }
        else if (option == "--tpc-place-x-cm")
        {
          options.map.tpc_place_x_cm =
              double_after(index, option);
        }
        else if (option == "--tpc-place-y-cm")
        {
          options.map.tpc_place_y_cm =
              double_after(index, option);
        }
        else if (option == "--tpc-place-z-cm")
        {
          options.map.tpc_place_z_cm =
              double_after(index, option);
        }
        else if (option == "--tpc-rotation-x-rad")
        {
          options.map.tpc_rotation_x_rad =
              double_after(index, option);
        }
        else if (option == "--tpc-rotation-y-rad")
        {
          options.map.tpc_rotation_y_rad =
              double_after(index, option);
        }
        else if (option == "--tpc-rotation-z-rad")
        {
          options.map.tpc_rotation_z_rad =
              double_after(index, option);
        }
        else if (option == "--global-tag")
        {
          options.global_tag = value_after(index, option);
        }
        else if (option == "--timestamp")
        {
          options.timestamp = uint64_after(index, option);
        }
        else if (option == "--label")
        {
          options.map.configuration_label = value_after(index, option);
        }
        else if (option == "--force")
        {
          options.map.overwrite = true;
        }
        else if (option == "--verbosity")
        {
          options.verbosity = int_after(index, option);
        }
        else
        {
          throw std::runtime_error("Unknown option " + option);
        }
      }
    }
    catch (const std::exception& exception)
    {
      error << exception.what() << "\n";
      return false;
    }

    if (options.help_requested)
    {
      return true;
    }
    if (options.adaptive)
    {
      if (options.map.interpolation_tolerance_um == 0.)
      {
        options.map.interpolation_tolerance_um = options.adaptive_policy.guard_um;
      }
      options.adaptive_policy.guard_um = options.map.interpolation_tolerance_um;
      if (!options.adaptive_end)
      {
        options.adaptive_end = PHGarfieldAdaptiveMap::RootCount(PHGarfieldAdaptiveMap::FromConfig(options.map));
      }
      if (options.adaptive_first >= options.adaptive_end || options.adaptive_end > PHGarfieldAdaptiveMap::RootCount(PHGarfieldAdaptiveMap::FromConfig(options.map)) ||
          options.map.phi_chunk_index != 0 || options.map.phi_chunk_count != 1 || !options.map.probe_cell_centers)
      {
        error << "Adaptive mode requires a nonempty valid root range, no phi chunks and the complete stencil\n";
        return false;
      }
    }
    else if (options.adaptive_range_set)
    {
      error << "Adaptive root ranges require --adaptive\n";
      return false;
    }
    if (options.map.output_filename.empty())
    {
      error << "--output is required\n";
      return false;
    }
    if (options.electric_field_map.empty() && options.electric_field_map_3d[0].empty() &&
        !options.allow_no_electric_corrections)
    {
      error << "--electric-field-map or both --electric-field-map-3d-side* are required unless "
               "--allow-no-electric-corrections is used\n";
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
  if (!require_regular_file(
          options.electric_field_map,
          "--electric-field-map") ||
      !require_regular_file(
          options.magnetic_field_map,
          "--magnetic-field-map"))
  {
    return 2;
  }
  for (const auto& paths : {options.electric_field_map_3d, options.frame_field_map_3d})
  {
    if (paths[0].empty() != paths[1].empty())
    {
      std::cerr << "Side-specific 3D payloads must be pinned for both sides\n";
      return 2;
    }
    for (const auto& path : paths)
    {
      if (!require_regular_file(path, "3D field payload"))
      {
        return 2;
      }
    }
  }
  if (!options.electric_field_map.empty() && !options.electric_field_map_3d[0].empty())
  {
    std::cerr << "Select either axisymmetric 2D or per-side 3D space-charge maps\n";
    return 2;
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
  if (options.electric_field_map.empty() && options.electric_field_map_3d[0].empty() &&
      (options.space_charge_scale_side0 != 0.0 ||
       options.space_charge_scale_side1 != 0.0))
  {
    std::cerr
        << "A nonzero --scale-side0 or --scale-side1 requires a valid "
           "--electric-field-map. For an explicit nominal E/B-only map, use "
           "--allow-no-electric-corrections together with zero scales on "
           "both sides.\n";
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

  auto garfield = std::make_unique<PHGarfield>(
      "PHGarfieldForwardMapGenerator",
      options.electric_field_map,
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
  if (!options.electric_field_map.empty() &&
      !garfield->HasElectricFieldCorrections())
  {
    std::cerr << "The requested electric-field correction map was not loaded\n";
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
  if (options.adaptive)
  {
    PHGarfieldAdaptiveMap::BuildSummary summary;
    return PHGarfieldAdaptiveMap::Generate(*garfield, options.map, options.adaptive_policy,
                                           options.adaptive_first, options.adaptive_end, summary, std::cout, std::cerr)
               ? 0
               : 1;
  }
  PHGarfieldForwardMap::Summary summary;
  if (!PHGarfieldForwardMap::Generate(
          *garfield,
          options.map,
          summary,
          std::cout,
          std::cerr))
  {
    return 1;
  }
  return 0;
}
catch (const std::exception& error)
{
  std::cerr << "Forward-map generation failed: " << error.what() << '\n';
  return 1;
}

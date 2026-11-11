#include "PHGarfieldMapTools.h"
#include "PHGarfieldRootOutput.h"

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TParameter.h>

#include <unistd.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PHGarfieldMapTools
{
  const std::string MetadataPath = "GarfieldFullTransport/Metadata";
  const std::array<const char*, 2> ElectricFieldMapPayloads = {{"Tpc_PolySeeding_EField3D_Side0", "Tpc_PolySeeding_EField3D_Side1"}};
  const std::array<const char*, 2> FrameFieldMapPayloads = {{"Tpc_PolySeeding_ModuleFrames3D_Side0", "Tpc_PolySeeding_ModuleFrames3D_Side1"}};

  void require(bool value, const std::string& message)
  {
    if (!value)
    {
      throw std::runtime_error(message);
    }
  }

  PHGarfieldCompactMap::Limits generation_limits(const Map::Policy& policy)
  {
    using Compact = PHGarfieldCompactMap;
    Compact::Limits limits;
    limits.max_nodes = std::min(policy.max_total_cells, Compact::HardMaximumNodes);
    limits.max_corners = std::min(policy.max_total_cells, Compact::HardMaximumCorners);
    limits.max_values = std::min(policy.max_total_points, Compact::HardMaximumValues);
    limits.max_bytes = std::min(policy.max_estimated_bytes, Compact::HardMaximumBytes);
    return limits;
  }

  std::array<std::string, 2> resolve_field_maps(
      const std::array<std::string, 2>& overrides, const std::array<const char*, 2>& payloads,
      const std::function<std::string(const std::string&)>& lookup)
  {
    auto resolved = overrides;
    for (std::size_t side = 0; side < resolved.size(); ++side)
    {
      if (resolved[side].empty())
      {
        resolved[side] = lookup(payloads[side]);
        require(!resolved[side].empty(), std::string("missing CDB payload ") + payloads[side]);
      }
    }
    return resolved;
  }

  Options options(int argc, char** argv, std::initializer_list<std::string> allowed)
  {
    Options result;
    for (int i = 1; i < argc; i += 2)
    {
      require(i + 1 < argc && result.emplace(argv[i], argv[i + 1]).second, "options require unique key/value pairs");
      require(std::find(allowed.begin(), allowed.end(), argv[i]) != allowed.end(), std::string("unsupported option ") + argv[i]);
    }
    return result;
  }

  std::uint64_t count(TDirectory& directory, const char* name)
  {
    const auto value = parameter<Long64_t>(directory, name);
    require(value >= 0, std::string("negative metadata ") + name);
    return static_cast<std::uint64_t>(value);
  }

  std::string named(TDirectory& directory, const char* name)
  {
    const auto* object = dynamic_cast<TNamed*>(directory.Get(name));
    require(object != nullptr && std::string(object->ClassName()) == "TNamed",
            std::string("missing/incorrect named metadata ") + name);
    return object->GetTitle();
  }

  void publish_exclusive(const std::string& temporary, const std::string& output)
  {
    require(::link(temporary.c_str(), output.c_str()) == 0, "exclusive publication failed; temporary evidence preserved");
    require(::unlink(temporary.c_str()) == 0, "output published but temporary unlink failed");
  }

  void write_compact_output(const std::string& path, const PHGarfieldCompactMap::Data& data,
                            const PHGarfieldCompactMap::Limits& limits, std::ostream& errors,
                            const std::function<void(TDirectory&)>& write_metadata)
  {
    using Compact = PHGarfieldCompactMap;
    TFile file(path.c_str(), "CREATE");
    require(!file.IsZombie(), "cannot create unique temporary ROOT output");
    auto* top = file.mkdir("GarfieldFullTransport");
    require(top, "cannot create transport directory");
    auto* metadata = top->mkdir("Metadata");
    require(metadata, "cannot create transport metadata");
    {
      TDirectory::TContext context(metadata);
      write_metadata(*metadata);
    }
    require(Compact::Write(top->mkdir("Runtime"), data, limits, errors), "compact array write failed");
    require(PHGarfieldFlushAndCloseRootOutput(file, errors), "ROOT flush/close failed");
  }

  void readback_compact_output(const std::string& path, PHGarfieldCompactMap::Data& data,
                               const PHGarfieldCompactMap::Limits& limits, std::ostream& errors,
                               bool complete)
  {
    using Compact = PHGarfieldCompactMap;
    const auto geometry_value = data.geometry;
    const auto content = Compact::ContentFingerprint(data);
    data = {};
    auto file = open(path);
    require(Compact::Read(file->GetDirectory("GarfieldFullTransport/Runtime"), geometry_value, data, limits, errors, complete),
            "independent ROOT readback failed");
    require(Compact::ContentFingerprint(data) == content, "ROOT readback changed content");
  }
  Map::Geometry geometry(TDirectory& directory)
  {
    Map::Geometry result;
    const auto positive = [&](const char* name)
    {
      const int value = parameter<int>(directory, name);
      require(value > 0, std::string("nonpositive geometry ") + name);
      return static_cast<unsigned int>(value);
    };
    result.phi_nodes = positive("phi_nodes");
    result.radial_nodes = positive("radial_nodes");
    result.longitudinal_nodes = positive("longitudinal_nodes");
    result.inner_radius_cm = parameter<double>(directory, "gas_inner_radius_cm");
    result.outer_radius_cm = parameter<double>(directory, "gas_outer_radius_cm");
    result.central_membrane_half_width_cm = parameter<double>(directory, "central_membrane_half_width_cm");
    result.readout_abs_z_cm = parameter<double>(directory, "readout_abs_z_cm");
    result.nominal_velocity_cm_ns = parameter<double>(directory, "nominal_phg4_drift_velocity_cm_per_ns");
    return result;
  }

  std::unique_ptr<TFile> open(const std::string& path)
  {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    require(file && !file->IsZombie() && !file->TestBit(TFile::kRecovered), "unreadable/recovered shard " + path);
    return file;
  }

  const std::array<const char*, 5> GenerationCounters = {{"adaptive_sampler_queries", "adaptive_reused_points", "adaptive_certified_leaves",
                                                          "adaptive_direct_leaves", "adaptive_budget_leaves"}};

  template <class T>
  bool parameter_signature(const TObject& object, std::ostringstream& output)
  {
    if (const auto* value = dynamic_cast<const TParameter<T>*>(&object))
    {
      output << value->GetVal();
      return true;
    }
    return false;
  }

  std::string signature(const TObject& object)
  {
    std::ostringstream output;
    output << object.ClassName() << ':' << object.GetName() << ':' << std::hexfloat << std::setprecision(17);
    if (parameter_signature<double>(object, output) || parameter_signature<int>(object, output) ||
        parameter_signature<Long64_t>(object, output))
    {
      return output.str();
    }
    const auto* text = dynamic_cast<const TNamed*>(&object);
    require(text && std::string(object.ClassName()) == "TNamed", "unsupported common metadata type " + std::string(object.ClassName()));
    output << text->GetTitle();
    return output.str();
  }

  // These fields describe execution/aggregation, not the common physical model.
  const std::set<std::string> VariableMetadata = {
      "generation_root_first", "generation_root_end", "generation_shards",
      "maximum_drift_time_ns", "reach_probe_trajectories", "adaptive_sampler_queries",
      "adaptive_reused_points", "adaptive_certified_leaves", "adaptive_direct_leaves", "adaptive_budget_leaves"};

  std::map<std::string, std::string> common_metadata(TDirectory& directory, Objects* saved)
  {
    std::map<std::string, std::string> result;
    std::set<std::string> seen;
    for (auto* entry : *directory.GetListOfKeys())
    {
      auto* key = dynamic_cast<TKey*>(entry);
      require(key != nullptr, "non-key in metadata list");
      const std::string name = key->GetName();
      require(seen.insert(name).second, "duplicate metadata key/cycle " + name);
      std::unique_ptr<TObject> object(key->ReadObj());
      require(object != nullptr, "unreadable metadata " + name);
      if (VariableMetadata.contains(name))
      {
        continue;
      }
      result[name] = signature(*object);
      if (saved)
      {
        saved->emplace(name, std::move(object));
      }
    }
    require(!result.empty(), "empty common metadata");
    return result;
  }

}  // namespace PHGarfieldMapTools

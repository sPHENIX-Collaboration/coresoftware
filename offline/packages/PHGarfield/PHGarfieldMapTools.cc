#include "PHGarfieldMapTools.h"
#include <cerrno>
#include "PHGarfieldCompactMap.h"
#include "PHGarfieldMapFingerprint.h"

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TParameter.h>

#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PHGarfieldMapTools
{
  const std::string MetadataPath = "GarfieldFullTransport/Metadata";
  const std::string AdaptivePath = "GarfieldFullTransport/Adaptive";

  void require(bool value, const std::string& message)
  {
    if (!value)
    {
      throw std::runtime_error(message);
    }
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

  std::string fingerprint(const std::string& text)
  {
    return PHGarfieldForwardMap::Fingerprint(text);
  }

  // Explicit output resources do not change any generation/refinement policy.
  std::uint64_t integer(const std::string& text)
  {
    require(!text.empty() && std::all_of(text.begin(), text.end(), [](char c)
                                         { return c >= '0' && c <= '9'; }),
            "nonnegative decimal count required");
    std::size_t consumed = 0;
    const auto value = std::stoull(text, &consumed);
    require(consumed == text.size() && value <= static_cast<std::uint64_t>(std::numeric_limits<Long64_t>::max()), "count exceeds signed ROOT metadata bound");
    return value;
  }
  void apply_capacity(Map::Policy& policy, const Capacity& capacity)
  {
    if (!capacity.enabled)
    {
      return;
    }
    require(capacity.points >= policy.max_total_points && capacity.cells >= policy.max_total_cells && capacity.bytes >= policy.max_estimated_bytes,
            "output capacity must monotonically preserve all original input resource limits");
    policy.max_total_points = capacity.points;
    policy.max_total_cells = capacity.cells;
    policy.max_estimated_bytes = capacity.bytes;
  }
  std::string capacity_prefix(const std::string& canonical)
  {
    const auto split = canonical.rfind(";max_points=");
    require(split != std::string::npos, "missing canonical adaptive capacity fields");
    return canonical.substr(0, split);
  }
  std::array<std::uint64_t, 3> canonical_capacity(const std::string& canonical)
  {
    const std::array<std::string, 3> names{{";max_points=", ";max_cells=", ";max_bytes="}};
    std::array<std::uint64_t, 3> values{};
    auto cursor = canonical.rfind(names[0]);
    require(cursor != std::string::npos, "missing canonical adaptive capacity");
    for (unsigned int index = 0; index < names.size(); ++index)
    {
      require(canonical.compare(cursor, names[index].size(), names[index]) == 0, "malformed canonical adaptive capacity");
      const auto first = cursor + names[index].size();
      const auto end = canonical.find(';', first);
      values[index] = integer(canonical.substr(first, end - first));
      require(values[index] > 0, "zero canonical capacity");
      if (index + 1 == names.size())
      {
        require(end == std::string::npos, "extra canonical capacity suffix");
      }
      else
      {
        require(end != std::string::npos, "truncated canonical capacity");
        cursor = end;
      }
    }
    return values;
  }
  // Same primitive byte order as the native content digest, excluding only its
  // three allocation ceilings. This is linear and never builds a giant string.
  std::string primitive_fingerprint(const Map::Data& data)
  {
    return PHGarfieldMapFingerprint::Adaptive(data, capacity_prefix(Map::BuildCanonicalPolicy(data)));
  }
  bool same_values(const std::array<double, 3>& first, const std::array<double, 3>& second)
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (std::bit_cast<std::uint64_t>(first[i]) != std::bit_cast<std::uint64_t>(second[i]))
      {
        return false;
      }
    }
    return true;
  }
  std::string quote(const std::string& value)
  {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value)
    {
      if (c == '"' || c == '\\')
      {
        out << '\\' << static_cast<char>(c);
      }
      else if (c < 32)
      {
        out << "\\u" << std::hex << std::setfill('0') << std::setw(4) << static_cast<unsigned int>(c) << std::dec;
      }
      else
      {
        out << static_cast<char>(c);
      }
    }
    out << '"';
    return out.str();
  }
  void write_exclusive(const std::string& path, const std::string& text)
  {
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
    require(fd >= 0, "cannot create exclusive proof temporary");
    std::size_t done = 0;
    while (done < text.size())
    {
      const auto count = ::write(fd, text.data() + done, text.size() - done);
      if (count < 0 && errno == EINTR)
      {
        continue;
      }
      if (count <= 0)
      {
        ::close(fd);
        throw std::runtime_error("proof write failed; partial evidence preserved");
      }
      done += static_cast<std::size_t>(count);
    }
    const bool synced = ::fsync(fd) == 0;
    const bool closed = ::close(fd) == 0;
    require(synced && closed, "proof flush/close failed");
  }
  void check_capacity_envelope(TDirectory& metadata, const Map::Data& data)
  {
    const std::array<const char*, 5> fields{{"merge_capacity_transform", "merge_input_adaptive_policy", "merge_output_adaptive_policy",
                                             "merge_capacity_independent_fingerprint_fnv1a64", "merge_source_shards_sha256_tsv"}};
    unsigned int present = 0;
    for (const auto* name : fields)
    {
      present += metadata.Get(name) != nullptr;
    }
    if (!present)
    {
      return;
    }
    require(present == fields.size() && named(metadata, fields[0]) == "monotone_aggregate_capacity_only_v1", "incomplete/invalid capacity transform provenance");
    const auto input = named(metadata, fields[1]);
    const auto output = named(metadata, fields[2]);
    require(output == Map::BuildCanonicalPolicy(data) && capacity_prefix(input) == capacity_prefix(output), "capacity transform changed nonresource adaptive policy");
    const auto before = canonical_capacity(input);
    const auto after = canonical_capacity(output);
    for (unsigned int i = 0; i < 3; ++i)
    {
      require(after[i] >= before[i], "capacity provenance is not monotone");
    }
    require(named(metadata, fields[3]) == primitive_fingerprint(data), "capacity-independent primitive fingerprint mismatch");
    require(!named(metadata, fields[4]).empty(), "capacity transform lacks original source SHA inventory");
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

  void check_summary(TDirectory& metadata, const Map::Data& data)
  {
    require(count(metadata, "generation_root_first") == data.first_root &&
                count(metadata, "generation_root_end") == data.end_root,
            "top metadata/adaptive root range mismatch");
    require(parameter<double>(metadata, "interpolation_tolerance_um") == data.policy.guard_um &&
                parameter<double>(metadata, "maximum_drift_time_ns") == data.maximum_drift_time_ns,
            "top metadata/adaptive guard or duration mismatch");
  }

  const std::array<const char*, 5> GenerationCounters = {{"adaptive_sampler_queries", "adaptive_reused_points", "adaptive_certified_leaves",
                                                          "adaptive_direct_leaves", "adaptive_budget_leaves"}};

  std::array<std::uint64_t, 5> generation_counters(TDirectory& metadata, const Map::Data& data)
  {
    std::array<std::uint64_t, 5> result{};
    for (unsigned int i = 0; i < result.size(); ++i)
    {
      result[i] = count(metadata, GenerationCounters[i]);
    }
    require(named(metadata, "reach_probe_trajectories") == std::to_string(result[0]),
            "adaptive trajectory/sampler counter mismatch");
    std::uint64_t certified = 0;
    std::uint64_t direct = 0;
    std::uint64_t budget = 0;
    for (const auto& cell : data.cells)
    {
      if (cell.split_axis >= 0)
      {
        continue;
      }
      if (cell.leaf_class == static_cast<int>(Map::LeafClass::Certified))
      {
        ++certified;
      }
      else
      {
        ++direct;
      }
      if (cell.fallback_reason >= static_cast<int>(Map::FallbackReason::DepthBudget) &&
          cell.fallback_reason <= static_cast<int>(Map::FallbackReason::QueryBudget))
      {
        ++budget;
      }
    }
    require(result[2] == certified && result[3] == direct && result[4] == budget,
            "adaptive generation leaf counters disagree with tree");
    return result;
  }

  void check_file(const std::string& path, bool complete)
  {
    auto file = open(path);
    auto* metadata = file->GetDirectory(MetadataPath.c_str());
    require(metadata && parameter<int>(*metadata, "schema_version") == Map::SchemaVersion,
            "missing schema4 metadata");
    require(named(*metadata, "schema_name") == (complete ? PHGarfieldForwardMap::SchemaName : "PHG4TpcGarfieldFullTransportLUTShard"), "wrong full/shard identity");
    const auto canonical = named(*metadata, "canonical_configuration");
    require(!canonical.empty() && fingerprint(canonical) == named(*metadata, "configuration_fingerprint_fnv1a64"),
            "common configuration fingerprint mismatch");
    Map::Data data;
    require(Map::Read(file->GetDirectory(AdaptivePath.c_str()), geometry(*metadata), data, std::cerr, complete),
            "adaptive content/structure failed validation");
    check_summary(*metadata, data);
    check_capacity_envelope(*metadata, data);
    (void) generation_counters(*metadata, data);
    std::uint64_t certified = 0;
    std::uint64_t direct = 0;
    for (const auto& cell : data.cells)
    {
      if (cell.split_axis >= 0)
      {
        continue;
      }
      if (cell.leaf_class == static_cast<int>(Map::LeafClass::Certified))
      {
        ++certified;
      }
      else
      {
        ++direct;
      }
    }
    std::cout << "{\"status\":\"PASS_STRUCTURE_ONLY\",\"first_root\":" << data.first_root
              << ",\"end_root\":" << data.end_root << ",\"points\":" << data.points.size()
              << ",\"cells\":" << data.cells.size() << ",\"roots\":" << data.roots.size()
              << ",\"certified_leaves\":" << certified << ",\"direct_leaves\":" << direct
              << ",\"guard_um\":" << data.policy.guard_um << ",\"content_fingerprint\":\""
              << Map::ContentFingerprint(data) << "\"}\n";
  }

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
      "schema_name", "generation_phi_chunk_index", "generation_phi_chunk_count",
      "generation_root_first", "generation_root_end", "adaptive_root_first", "adaptive_root_end",
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

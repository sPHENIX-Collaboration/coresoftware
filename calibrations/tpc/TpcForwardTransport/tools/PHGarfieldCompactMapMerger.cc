#include "PHGarfieldCompactMap.h"
#include "PHGarfieldMapTools.h"

#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>

#include <unistd.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
  using PHGarfieldMapTools::common_metadata;
  using PHGarfieldMapTools::count;
  using PHGarfieldMapTools::GenerationCounters;
  using PHGarfieldMapTools::geometry;
  using PHGarfieldMapTools::MetadataPath;
  using PHGarfieldMapTools::named;
  using PHGarfieldMapTools::Objects;
  using PHGarfieldMapTools::open;
  using PHGarfieldMapTools::options;
  using PHGarfieldMapTools::parameter;
  using PHGarfieldMapTools::publish_exclusive;
  using PHGarfieldMapTools::readback_compact_output;
  using PHGarfieldMapTools::require;
  using PHGarfieldMapTools::write_compact_output;
  using Compact = PHGarfieldCompactMap;

  struct Input
  {
    std::string path;
    std::uintmax_t bytes;
    std::filesystem::file_time_type modified;
  };

  void unchanged(const Input& input)
  {
    require(std::filesystem::file_size(input.path) == input.bytes &&
                std::filesystem::last_write_time(input.path) == input.modified,
            "input changed during merge: " + input.path);
  }

  void add(std::uint64_t& total, std::uint64_t amount, std::uint64_t maximum)
  {
    require(total <= maximum && amount <= maximum - total, "collection count exceeds limit");
    total += amount;
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--help")
  {
    std::cout << "Usage: PHGarfieldCompactMapMerger --input-list files.txt --output runtime.root [--complete yes|no]\n"
                 "List one LUT path per line, in ascending contiguous root-range order.\n"
                 "The default requires complete coverage; --complete no writes a partial collection for a later merge.\n";
    return 0;
  }
  try
  {
    const auto arguments = options(argc, argv, {"--input-list", "--output", "--complete"});
    const auto get = [&](const char* key) -> const std::string&
    {
      require(arguments.contains(key), std::string("missing ") + key);
      return arguments.at(key);
    };
    const auto mode = arguments.find("--complete");
    require(mode == arguments.end() || mode->second == "yes" || mode->second == "no", "complete must be yes or no");
    const bool complete = mode == arguments.end() || mode->second == "yes";
    const auto& output = get("--output");
    require(!std::filesystem::exists(output), "output already exists");
    const auto list_path = std::filesystem::absolute(get("--input-list"));
    std::ifstream list(list_path);
    require(list.is_open(), "input list unreadable");
    Compact::Limits limits;
    std::uint64_t total_nodes = 0;
    std::uint64_t total_corners = 0;
    std::uint64_t total_values = 0;
    std::uint64_t total_roots = 0;
    std::vector<Input> inputs;
    std::set<std::string> unique;
    std::string line;
    while (std::getline(list, line))
    {
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      if (line.empty())
      {
        continue;
      }
      const std::filesystem::path entry(line);
      const auto path = std::filesystem::canonical(entry.is_absolute() ? entry : list_path.parent_path() / entry).string();
      require(unique.insert(path).second && inputs.size() < 1000000, "duplicate input or excessive input count");
      Input input{path, std::filesystem::file_size(path), std::filesystem::last_write_time(path)};
      auto file = open(path);
      auto* runtime = file->GetDirectory("GarfieldFullTransport/Runtime");
      require(runtime, "missing runtime arrays");
      add(total_nodes, count(*runtime, "node_count"), limits.max_nodes);
      add(total_corners, count(*runtime, "corner_record_count"), limits.max_corners);
      add(total_values, count(*runtime, "value_count"), limits.max_values);
      total_roots = std::max(total_roots, count(*runtime, "root_lookup_count"));
      unchanged(input);
      inputs.push_back(input);
    }
    require(list.eof() && !inputs.empty(), "empty or unreadable input list");
    std::uint64_t bytes = 0;
    require(total_roots <= limits.max_roots &&
                Compact::EstimatedBytes(total_nodes, total_corners, total_values, total_roots, bytes) && bytes <= limits.max_bytes,
            "collection allocation budget exceeded; merge smaller groups first");
    Compact::Data data;
    data.nodes.reserve(total_nodes);
    data.corners.reserve(total_corners);
    data.values.reserve(total_values);
    Objects saved;
    std::map<std::string, std::string> expected;
    std::array<std::uint64_t, 5> counters{};
    std::uint64_t shard_count = 0;
    for (const auto& input : inputs)
    {
      unchanged(input);
      auto file = open(input.path);
      auto* metadata = file->GetDirectory(MetadataPath.c_str());
      require(metadata && parameter<int>(*metadata, "schema_version") == Compact::SchemaVersion &&
                  named(*metadata, "schema_name") == Compact::SchemaName,
              "missing or unsupported LUT metadata");
      const auto common = common_metadata(*metadata, expected.empty() ? &saved : nullptr);
      require(expected.empty() || expected == common, "physical/refinement metadata mismatch");
      Compact::Data part;
      require(Compact::Read(file->GetDirectory("GarfieldFullTransport/Runtime"), geometry(*metadata), part, limits, std::cerr, false),
              "invalid LUT shard");
      require(PHGarfieldForwardMap::Fingerprint(named(*metadata, "canonical_configuration")) == part.configuration_fingerprint_fnv1a64 &&
                  named(*metadata, "configuration_fingerprint_fnv1a64") == part.configuration_fingerprint_fnv1a64,
              "physical configuration fingerprint mismatch");
      require(count(*metadata, "generation_root_first") == part.first_root && count(*metadata, "generation_root_end") == part.end_root &&
                  parameter<double>(*metadata, "maximum_drift_time_ns") == part.maximum_drift_time_ns &&
                  parameter<double>(*metadata, "interpolation_tolerance_um") == part.guard_um,
              "metadata/runtime range, duration or guard mismatch");
      if (expected.empty())
      {
        expected = common;
      }
      require(Compact::Append(data, part, limits, std::cerr), "incompatible, overlapping, missing or out-of-order root range");
      for (std::size_t index = 0; index < counters.size(); ++index)
      {
        add(counters[index], count(*metadata, GenerationCounters[index]), std::numeric_limits<Long64_t>::max());
      }
      add(shard_count, count(*metadata, "generation_shards"), std::numeric_limits<Long64_t>::max());
      unchanged(input);
    }
    require(data.nodes.size() == total_nodes && data.corners.size() == total_corners && data.values.size() == total_values,
            "preflight/actual array counts differ");
    require(Compact::Validate(data, limits, std::cerr, complete), "invalid collection or incomplete requested coverage");
    const auto first = data.first_root;
    const auto end = data.end_root;
    const auto temporary = output + ".partial." + std::to_string(getpid());
    write_compact_output(temporary, data, limits, std::cerr, [&](TDirectory&)
                         {
      for (const auto& [name, object] : saved)
      {
        (void) name;
        require(object->Write() > 0, "common metadata write failed");
      }
      TParameter<Long64_t>("generation_root_first", first).Write();
      TParameter<Long64_t>("generation_root_end", end).Write();
      TParameter<Long64_t>("generation_shards", shard_count).Write();
      TParameter<double>("maximum_drift_time_ns", data.maximum_drift_time_ns).Write();
      TNamed("reach_probe_trajectories", std::to_string(counters[0]).c_str()).Write();
      for (std::size_t index = 0; index < counters.size(); ++index)
      {
        TParameter<Long64_t>(GenerationCounters[index], counters[index]).Write();
      } });
    readback_compact_output(temporary, data, limits, std::cerr, complete);
    for (const auto& input : inputs)
    {
      unchanged(input);
    }
    publish_exclusive(temporary, output);
    std::cout << "Wrote adaptive LUT " << output << " roots=[" << first << ',' << end << ") from " << inputs.size()
              << " inputs; nodes=" << data.nodes.size() << " corners=" << data.corners.size() << " values=" << data.values.size() << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "LUT merge failed: " << error.what() << '\n';
    return 1;
  }
}

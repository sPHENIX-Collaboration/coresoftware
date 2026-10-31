#include "PHGarfieldCompactMap.h"
#include "PHGarfieldMapTools.h"

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

namespace
{
  using PHGarfieldMapTools::AdaptivePath;
  using PHGarfieldMapTools::apply_capacity;
  using PHGarfieldMapTools::Capacity;
  using PHGarfieldMapTools::check_capacity_envelope;
  using PHGarfieldMapTools::check_file;
  using PHGarfieldMapTools::check_summary;
  using PHGarfieldMapTools::common_metadata;
  using PHGarfieldMapTools::count;
  using PHGarfieldMapTools::fingerprint;
  using PHGarfieldMapTools::generation_counters;
  using PHGarfieldMapTools::GenerationCounters;
  using PHGarfieldMapTools::geometry;
  using PHGarfieldMapTools::integer;
  using PHGarfieldMapTools::Map;
  using PHGarfieldMapTools::MetadataPath;
  using PHGarfieldMapTools::named;
  using PHGarfieldMapTools::Objects;
  using PHGarfieldMapTools::open;
  using PHGarfieldMapTools::parameter;
  using PHGarfieldMapTools::primitive_fingerprint;
  using PHGarfieldMapTools::quote;
  using PHGarfieldMapTools::require;
  using PHGarfieldMapTools::same_values;
  using PHGarfieldMapTools::write_exclusive;

  void append(Map::Data& merged, const Map::Data& shard,
              std::map<std::array<double, 3>, std::uint64_t>& coordinates, const Capacity& capacity)
  {
    require(shard.first_root == merged.end_root && shard.end_root > shard.first_root,
            "missing/overlapping/empty adaptive root range");
    require(Map::BuildCanonicalPolicy(merged) == Map::BuildCanonicalPolicy(shard), "adaptive geometry/policy mismatch");
    auto output_policy = merged.policy;
    apply_capacity(output_policy, capacity);
    std::uint64_t estimated = 0;
    // Check a conservative upper bound before allocating, even when point reuse
    // would make the final result smaller. Failures retain their evidence.
    require(shard.points.size() <= output_policy.max_total_points - merged.points.size() &&
                shard.cells.size() <= output_policy.max_total_cells - merged.cells.size() &&
                Map::EstimatedBytes(merged.points.size() + shard.points.size(), merged.cells.size() + shard.cells.size(),
                                    Map::RootCount(merged.geometry), estimated) &&
                estimated <= output_policy.max_estimated_bytes,
            "merged adaptive allocation budget exceeded");
    std::vector<std::uint64_t> point_ids;
    point_ids.reserve(shard.points.size());
    for (auto point : shard.points)
    {
      const auto found = coordinates.find(point.coordinate);
      if (found == coordinates.end())
      {
        point.id = merged.points.size();
        point_ids.push_back(point.id);
        coordinates.emplace(point.coordinate, point.id);
        merged.points.push_back(point);
      }
      else
      {
        const auto& previous = merged.points[found->second];
        require(previous.reaches == point.reaches && previous.forward_status == point.forward_status && same_values(previous.value, point.value),
                "duplicate coordinate has inconsistent physical transport");
        point_ids.push_back(found->second);
      }
    }
    const std::uint64_t offset = merged.cells.size();
    for (auto cell : shard.cells)
    {
      cell.id += offset;
      for (auto& corner : cell.corner_id)
      {
        if (corner >= 0)
        {
          corner = point_ids.at(corner);
        }
      }
      for (auto& child : cell.child_id)
      {
        if (child >= 0)
        {
          child += offset;
        }
      }
      merged.cells.push_back(cell);
    }
    for (auto root : shard.roots)
    {
      root.cell_id += offset;
      merged.roots.push_back(root);
    }
    merged.end_root = shard.end_root;
    merged.maximum_drift_time_ns = std::max(merged.maximum_drift_time_ns, shard.maximum_drift_time_ns);
  }

  void merge(const std::string& destination, const std::vector<std::string>& paths, const Capacity& capacity)
  {
    require(!std::filesystem::exists(destination) && !paths.empty(), "output exists or no input shards");
    if (capacity.enabled)
    {
      require(!capacity.proof_output.empty() && !std::filesystem::exists(capacity.proof_output) &&
                  std::filesystem::absolute(capacity.proof_output).lexically_normal() != std::filesystem::absolute(destination).lexically_normal(),
              "capacity proof already exists or overlaps output");
      for (const auto& path : paths)
      {
        require(std::filesystem::absolute(capacity.proof_output).lexically_normal() != std::filesystem::absolute(path).lexically_normal(), "capacity proof overlaps source");
      }
    }
    struct Input
    {
      std::uint64_t first;
      std::string path;
    };
    std::vector<Input> inputs;
    std::set<std::string> unique_paths;
    for (const auto& path : paths)
    {
      require(unique_paths.insert(std::filesystem::canonical(path).string()).second, "duplicate input shard path");
      auto file = open(path);
      auto* directory = file->GetDirectory(AdaptivePath.c_str());
      require(directory != nullptr, "missing adaptive directory");
      inputs.push_back({count(*directory, "first_root"), path});
    }
    std::sort(inputs.begin(), inputs.end(), [](const Input& a, const Input& b)
              { return a.first < b.first; });
    Map::Data merged;
    Objects saved;
    std::map<std::string, std::string> expected;
    std::map<std::array<double, 3>, std::uint64_t> coordinates;
    std::string sources;
    std::string source_sha_inventory;
    std::ostringstream source_json;
    source_json << '[';
    bool first_source = true;
    std::array<std::uint64_t, 5> totals{};
    for (const auto& input : inputs)
    {
      const auto original_size = std::filesystem::file_size(input.path);
      const auto original_time = std::filesystem::last_write_time(input.path);
      const auto input_sha = capacity.enabled ? PHGarfieldCompactMap::FileSha256(input.path) : std::string{};
      auto file = open(input.path);
      auto* metadata = file->GetDirectory(MetadataPath.c_str());
      require(metadata && parameter<int>(*metadata, "schema_version") == Map::SchemaVersion &&
                  named(*metadata, "schema_name") == "PHG4TpcGarfieldFullTransportLUTShard",
              "not an adaptive LUT shard");
      const auto canonical = named(*metadata, "canonical_configuration");
      require(!canonical.empty() && fingerprint(canonical) == named(*metadata, "configuration_fingerprint_fnv1a64"),
              "common configuration fingerprint mismatch");
      Map::Data shard;
      require(Map::Read(file->GetDirectory(AdaptivePath.c_str()), geometry(*metadata), shard, std::cerr, false),
              "adaptive shard content/structure failed validation");
      check_summary(*metadata, shard);
      const auto counts = generation_counters(*metadata, shard);
      for (unsigned int i = 0; i < counts.size(); ++i)
      {
        require(counts[i] <= static_cast<std::uint64_t>(std::numeric_limits<Long64_t>::max()) - totals[i],
                "adaptive generation counter overflow");
        totals[i] += counts[i];
      }
      const auto common = common_metadata(*metadata, expected.empty() ? &saved : nullptr);
      if (expected.empty())
      {
        expected = common;
        merged.geometry = shard.geometry;
        merged.policy = shard.policy;
        auto policy_check = Map::Data{};
        policy_check.geometry = shard.geometry;
        policy_check.policy = shard.policy;
        apply_capacity(policy_check.policy, capacity);
        require(Map::Validate(policy_check, std::cerr, false), "invalid explicit output resource policy");
      }
      require(common == expected, "common physical/numerical metadata differs across shards");
      require(std::filesystem::file_size(input.path) == original_size && std::filesystem::last_write_time(input.path) == original_time,
              "source shard changed during validated read");
      append(merged, shard, coordinates, capacity);
      if (capacity.enabled)
      {
        const auto path = std::filesystem::canonical(input.path).string();
        const auto& sha = input_sha;
        source_sha_inventory += path + "\t" + sha + "\t" + Map::ContentFingerprint(shard) + "\n";
        if (!first_source)
        {
          source_json << ',';
        }
        first_source = false;
        source_json << "{\"path\":" << quote(path) << ",\"sha256\":" << quote(sha) << ",\"bytes\":" << std::filesystem::file_size(path)
                    << ",\"first_root\":" << shard.first_root << ",\"end_root\":" << shard.end_root << ",\"points\":" << shard.points.size()
                    << ",\"cells\":" << shard.cells.size() << ",\"content_fingerprint\":" << quote(Map::ContentFingerprint(shard)) << '}';
      }
      sources += std::filesystem::canonical(input.path).string() + "\t" + Map::ContentFingerprint(shard) + "\n";
    }
    source_json << ']';
    const auto input_policy = Map::BuildCanonicalPolicy(merged);
    const auto primitive_before = capacity.enabled ? primitive_fingerprint(merged) : std::string{};
    apply_capacity(merged.policy, capacity);
    require(!capacity.enabled || primitive_before == primitive_fingerprint(merged), "resource-only transform changed primitives or per-root policy");
    const auto output_policy = Map::BuildCanonicalPolicy(merged);
    require(Map::Validate(merged, std::cerr, true), "merged tree incomplete/invalid");
    coordinates.clear();
    const auto unique = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::string partial = destination + ".partial." + std::to_string(unique);
    {
      TFile output(partial.c_str(), "CREATE");
      require(!output.IsZombie(), "cannot create partial merged ROOT file");
      auto* transport = output.mkdir("GarfieldFullTransport");
      auto* metadata = transport->mkdir("Metadata");
      metadata->cd();
      for (const auto& item : saved)
      {
        require(item.second->Write() > 0, "metadata write failed");
      }
      TNamed("schema_name", PHGarfieldForwardMap::SchemaName).Write();
      TParameter<int>("generation_phi_chunk_index", 0).Write();
      TParameter<int>("generation_phi_chunk_count", 1).Write();
      TParameter<Long64_t>("generation_root_first", 0).Write();
      TParameter<Long64_t>("generation_root_end", Map::RootCount(merged.geometry)).Write();
      TParameter<double>("maximum_drift_time_ns", merged.maximum_drift_time_ns).Write();
      TNamed("reach_probe_trajectories", std::to_string(totals[0]).c_str()).Write();
      for (unsigned int i = 0; i < totals.size(); ++i)
      {
        TParameter<Long64_t>(GenerationCounters[i], static_cast<Long64_t>(totals[i])).Write();
      }
      TNamed("generation_adaptive_shards_tsv", sources.c_str()).Write();
      if (capacity.enabled)
      {
        TNamed("merge_capacity_transform", "monotone_aggregate_capacity_only_v1").Write();
        TNamed("merge_input_adaptive_policy", input_policy.c_str()).Write();
        TNamed("merge_output_adaptive_policy", output_policy.c_str()).Write();
        TNamed("merge_capacity_independent_fingerprint_fnv1a64", primitive_before.c_str()).Write();
        TNamed("merge_source_shards_sha256_tsv", source_sha_inventory.c_str()).Write();
      }
      auto* adaptive = transport->mkdir("Adaptive");
      require(Map::Write(adaptive, merged, std::cerr), "adaptive content write failed");
      output.Write();
      require(!output.TestBit(TFile::kWriteError), "ROOT write error");
      output.Close();
    }
    // Do not retain the merged data while independently reopening the complete
    // file. Runtime likewise owns one model, rather than all source ROOT files.
    const auto content = Map::ContentFingerprint(merged);
    const auto point_count = merged.points.size();
    const auto cell_count = merged.cells.size();
    const auto root_count = merged.roots.size();
    std::uint64_t actual_estimated_bytes = 0;
    require(Map::EstimatedBytes(point_count, cell_count, Map::RootCount(merged.geometry), actual_estimated_bytes), "merged estimated-byte arithmetic overflow");
    const auto expected_geometry = merged.geometry;
    merged = {};
    {
      auto reopened = open(partial);
      Map::Data checked;
      require(Map::Read(reopened->GetDirectory(AdaptivePath.c_str()), expected_geometry, checked, std::cerr, true) &&
                  Map::ContentFingerprint(checked) == content,
              "merged independent reopen failed");
      if (capacity.enabled)
      {
        require(primitive_fingerprint(checked) == primitive_before, "independent reopen changed capacity-independent primitives");
        auto* metadata = reopened->GetDirectory(MetadataPath.c_str());
        require(metadata && metadata->Get("merge_capacity_transform"), "missing reopened capacity metadata");
        check_capacity_envelope(*metadata, checked);
      }
    }
    std::filesystem::create_hard_link(partial, destination);
    std::filesystem::remove(partial);
    if (capacity.enabled)
    {
      std::ostringstream proof;
      proof << "{\"schema\":\"phgarfield_adaptive_capacity_merge_proof_v1\",\"status\":\"PASS_CAPACITY_ONLY_PRIMITIVES_PRESERVED\","
            << "\"transform\":\"monotone_aggregate_capacity_only_v1\",\"input_adaptive_policy\":" << quote(input_policy)
            << ",\"output_adaptive_policy\":" << quote(output_policy) << ",\"capacity_independent_fingerprint_fnv1a64\":" << quote(primitive_before)
            << ",\"content_fingerprint_fnv1a64\":" << quote(content)
            << ",\"output_capacity\":{\"points\":" << capacity.points << ",\"cells\":" << capacity.cells << ",\"estimated_bytes\":" << capacity.bytes << '}'
            << ",\"output_counts\":{\"points\":" << point_count << ",\"cells\":" << cell_count << ",\"roots\":" << root_count << ",\"estimated_bytes\":" << actual_estimated_bytes << '}'
            << ",\"sources\":" << source_json.str()
            << ",\"output\":{\"path\":" << quote(std::filesystem::canonical(destination).string()) << ",\"sha256\":" << quote(PHGarfieldCompactMap::FileSha256(destination))
            << ",\"bytes\":" << std::filesystem::file_size(destination) << "},\"numerical_values_recomputed\":false}\n";
      write_exclusive(capacity.proof_output, proof.str());
    }
    std::cout << "PASS: complete adaptive tree merged from " << paths.size() << " shards; content " << content
              << "; numerical/physics approval is separate\n";
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--help")
  {
    std::cout << "Usage: PHGarfieldAdaptiveMapMerger OUTPUT.root SHARD.root [SHARD.root ...]\n"
                 "       PHGarfieldAdaptiveMapMerger --check-shard FILE.root\n"
                 "       PHGarfieldAdaptiveMapMerger --check-complete FILE.root\n"
                 "Optional aggregate limits: --output-capacity-points N --output-capacity-cells N\n"
                 "  --output-capacity-bytes N --capacity-proof-output FILE.json\n";
    return 0;
  }
  try
  {
    if (argc == 3 && (std::string(argv[1]) == "--check-shard" || std::string(argv[1]) == "--check-complete"))
    {
      check_file(argv[2], std::string(argv[1]) == "--check-complete");
      return 0;
    }
    require(argc >= 3, "Usage: PHGarfieldAdaptiveMapMerger OUTPUT.root SHARD.root [SHARD.root ...]");
    Capacity capacity;
    std::map<std::string, std::string> options;
    int first = 1;
    while (first < argc && std::string(argv[first]).starts_with("--"))
    {
      require(first + 1 < argc && options.emplace(argv[first], argv[first + 1]).second, "capacity options require unique key/value pairs");
      first += 2;
    }
    const std::set<std::string> allowed{{"--output-capacity-points", "--output-capacity-cells", "--output-capacity-bytes", "--capacity-proof-output"}};
    for (const auto& [key, value] : options)
    {
      (void) value;
      require(allowed.contains(key), "unknown capacity option " + key);
    }
    if (!options.empty())
    {
      require(options.size() == 4, "capacity override requires all three output limits and a proof output");
      capacity.enabled = true;
      capacity.points = integer(options.at("--output-capacity-points"));
      capacity.cells = integer(options.at("--output-capacity-cells"));
      capacity.bytes = integer(options.at("--output-capacity-bytes"));
      require(capacity.points > 0 && capacity.cells > 0 && capacity.bytes >= 1024 && capacity.bytes <= Map::HardMaximumEstimatedBytes, "explicit capacity outside engineering limits");
      capacity.proof_output = options.at("--capacity-proof-output");
    }
    require(first + 1 < argc, "output and at least one source shard required");
    merge(argv[first], std::vector<std::string>(argv + first + 1, argv + argc), capacity);
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}

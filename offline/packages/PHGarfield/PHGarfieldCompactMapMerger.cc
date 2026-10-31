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

#include <TTree.h>
#include <fstream>
#include "PHGarfieldRootOutput.h"

namespace
{
  using PHGarfieldMapTools::AdaptivePath;
  using PHGarfieldMapTools::check_summary;
  using PHGarfieldMapTools::common_metadata;
  using PHGarfieldMapTools::count;
  using PHGarfieldMapTools::fingerprint;
  using PHGarfieldMapTools::GenerationCounters;
  using PHGarfieldMapTools::geometry;
  using PHGarfieldMapTools::Map;
  using PHGarfieldMapTools::MetadataPath;
  using PHGarfieldMapTools::named;
  using PHGarfieldMapTools::Objects;
  using PHGarfieldMapTools::open;
  using PHGarfieldMapTools::parameter;
  using PHGarfieldMapTools::quote;
  using PHGarfieldMapTools::require;
  using Compact = PHGarfieldCompactMap;
  struct Input
  {
    std::string path, sha, content;
    std::uintmax_t bytes;
    std::filesystem::file_time_type modified;
  };
  void unchanged(const Input& input)
  {
    require(std::filesystem::file_size(input.path) == input.bytes &&
                std::filesystem::last_write_time(input.path) == input.modified,
            "source changed during conversion");
  }
  void append_compact(Compact::Data& out, const Compact::Data& part)
  {
    require(part.first_root == out.end_root && part.end_root > part.first_root,
            "manifest shards must be in ascending contiguous first_root order (missing, duplicate, overlapping or out-of-order range)");
    require(part.source.configuration_fingerprint_fnv1a64 == out.source.configuration_fingerprint_fnv1a64 &&
                part.source.canonical_policy == out.source.canonical_policy && part.max_depth == out.max_depth &&
                part.max_nodes_per_root == out.max_nodes_per_root && part.max_queries_per_root == out.max_queries_per_root &&
                part.guard_um == out.guard_um,
            "incompatible shard physical/refinement policy");
    const auto no = out.nodes.size();
    const auto co = out.corners.size();
    const auto vo = out.values.size();
    const Compact::Limits limits;
    require(part.proof.source_points <= limits.max_source_points - out.proof.source_points, "source point count exceeds collection limit");
    require(no + part.nodes.size() < Compact::Invalid && co + part.corners.size() < Compact::Invalid &&
                vo + part.values.size() < Compact::Invalid,
            "compact offset overflow");
    for (auto n : part.nodes)
    {
      if ((n.tag & 3U) < 3U)
      {
        n.first += no;
        n.second += no;
      }
      else if (((n.tag >> 2U) & 7U) == static_cast<unsigned int>(Compact::LeafClass::Certified))
      {
        n.first += co;
      }
      out.nodes.push_back(n);
    }
    for (auto c : part.corners)
    {
      for (auto& v : c)
      {
        v += vo;
      }
      out.corners.push_back(c);
    }
    out.values.insert(out.values.end(), part.values.begin(), part.values.end());
    for (auto r = part.first_root; r < part.end_root; ++r)
    {
      require(out.root_lookup.at(r) == Compact::Invalid && part.root_lookup.at(r) != Compact::Invalid, "duplicate/missing compact root");
      out.root_lookup[r] = part.root_lookup[r] + no;
    }
    out.end_root = part.end_root;
    out.maximum_drift_time_ns = std::max(out.maximum_drift_time_ns, part.maximum_drift_time_ns);
    auto& p = out.proof;
    const auto& q = part.proof;
    p.source_points += q.source_points;
    p.source_cells += q.source_cells;
    p.roots_checked += q.roots_checked;
    p.bounds_checked += q.bounds_checked;
    p.routes_checked += q.routes_checked;
    p.certified_leaves_checked += q.certified_leaves_checked;
    p.corner_values_checked += q.corner_values_checked;
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--help")
  {
    std::cout << "Usage: PHGarfieldCompactMapMerger --mode convert|merge --input-list INPUTS.tsv\n"
                 "  --expect-input-sha256 SHA256 --output OUTPUT.root --summary OUTPUT.json\n"
                 "  --complete yes|no\n"
                 "INPUTS.tsv columns: absolute file path, SHA256, content fingerprint (tab-separated).\n";
    return 0;
  }
  try
  {
    std::map<std::string, std::string> options;
    for (int i = 1; i < argc; i += 2)
    {
      require(i + 1 < argc && options.emplace(argv[i], argv[i + 1]).second, "unique option/value pairs required");
    }
    const std::set<std::string> allowed{"--mode", "--input-list", "--expect-input-sha256", "--output", "--summary", "--complete"};
    for (const auto& [key, value] : options)
    {
      (void) value;
      require(allowed.contains(key), "unknown option " + key);
    }
    const auto get = [&](const char* key) -> const std::string&
    { require(options.contains(key),std::string("missing ")+key);return options.at(key); };
    const bool convert = get("--mode") == "convert";
    require(convert || get("--mode") == "merge", "mode must be convert or merge");
    require(get("--complete") == "yes" || get("--complete") == "no", "complete must be yes or no");
    const bool complete = get("--complete") == "yes";
    const std::string output = get("--output");
    const std::string summary = get("--summary");
    require(!std::filesystem::exists(output) && !std::filesystem::exists(summary) &&
                std::filesystem::absolute(output) != std::filesystem::absolute(summary),
            "preserve existing/distinct outputs");
    const auto executable_sha = Compact::FileSha256("/proc/self/exe");
    const auto input_sha = Compact::FileSha256(get("--input-list"));
    require(input_sha == get("--expect-input-sha256"), "input manifest checksum mismatch");
    std::ifstream list(get("--input-list"));
    std::ostringstream list_text;
    list_text << list.rdbuf();
    require(list.good(), "input manifest unreadable");
    const auto manifest = list_text.str();
    std::istringstream lines(manifest);
    std::string line;
    std::vector<Input> inputs;
    std::set<std::string> unique;
    Compact::Limits limits;
    std::uint64_t total_nodes = 0;
    std::uint64_t total_corners = 0;
    std::uint64_t total_values = 0;
    const auto add = [](std::uint64_t& total, std::uint64_t amount, std::uint64_t maximum)
    { require(total<=maximum && amount<=maximum-total,"declared collection count exceeds limit");total+=amount; };
    while (std::getline(lines, line))
    {
      std::istringstream fields(line);
      std::string path;
      std::string sha;
      std::string content;
      std::string extra;
      require(static_cast<bool>(std::getline(fields, path, '\t')) && static_cast<bool>(std::getline(fields, sha, '\t')) &&
                  static_cast<bool>(std::getline(fields, content, '\t')) && !std::getline(fields, extra, '\t'),
              "three TSV columns required");
      path = std::filesystem::canonical(path).string();
      require(unique.insert(path).second && inputs.size() < 1000000, "duplicate source or excessive input count");
      Input input{path, sha, content, std::filesystem::file_size(path), std::filesystem::last_write_time(path)};
      require(Compact::FileSha256(path) == sha, "source checksum mismatch: " + path);
      auto file = open(path);
      if (convert)
      {
        auto* metadata = file->GetDirectory(MetadataPath.c_str());
        require(metadata, "source metadata missing");
        auto* cells = dynamic_cast<TTree*>(file->Get((AdaptivePath + "/Cells").c_str()));
        auto* points = dynamic_cast<TTree*>(file->Get((AdaptivePath + "/Points").c_str()));
        require(cells && points && cells->GetEntries() >= 0 && points->GetEntries() >= 0, "missing source arrays");
        const auto certified = count(*metadata, "adaptive_certified_leaves");
        require(certified <= limits.max_corners, "excessive certified count");
        add(total_nodes, cells->GetEntries(), limits.max_nodes);
        add(total_corners, certified, limits.max_corners);
        add(total_values, std::min<std::uint64_t>(points->GetEntries(), 8 * certified), limits.max_values);
      }
      else
      {
        auto* runtime = file->GetDirectory("GarfieldFullTransport/Runtime");
        require(runtime, "missing compact arrays");
        add(total_nodes, count(*runtime, "node_count"), limits.max_nodes);
        add(total_corners, count(*runtime, "corner_record_count"), limits.max_corners);
        add(total_values, count(*runtime, "value_count"), limits.max_values);
      }
      unchanged(input);
      inputs.push_back(input);
    }
    require(!inputs.empty(), "empty input manifest");
    std::uint64_t estimate = 0;
    require(total_nodes <= limits.max_nodes && total_corners <= limits.max_corners && total_values <= limits.max_values &&
                Compact::EstimatedBytes(total_nodes, total_corners, total_values, limits.max_roots, estimate) && estimate <= limits.max_bytes,
            "collection allocation budget exceeded");
    Compact::Data data;
    data.nodes.reserve(total_nodes);
    data.corners.reserve(total_corners);
    data.values.reserve(total_values);
    Objects saved;
    std::map<std::string, std::string> expected;
    std::array<std::uint64_t, 5> counters{};
    for (const auto& input : inputs)
    {
      unchanged(input);
      auto file = open(input.path);
      auto* metadata = file->GetDirectory(convert ? MetadataPath.c_str() : "GarfieldFullTransport/SourceMetadata");
      require(metadata && parameter<int>(*metadata, "schema_version") == 4, "source metadata/version missing");
      const auto common = common_metadata(*metadata, expected.empty() ? &saved : nullptr);
      require(expected.empty() || expected == common, "physical/numerical metadata mismatch");
      Compact::Data part;
      if (convert)
      {
        require(named(*metadata, "schema_name") == "PHG4TpcGarfieldFullTransportLUTShard", "not a training shard");
        Map::Data source;
        require(Map::Read(file->GetDirectory(AdaptivePath.c_str()), geometry(*metadata), source, std::cerr, false), "invalid training shard");
        check_summary(*metadata, source);
        require(Map::ContentFingerprint(source) == input.content, "source content differs from audited manifest");
        Compact::SourceBinding binding{input.sha, named(*metadata, "configuration_fingerprint_fnv1a64"), input.content,
                                       Map::BuildCanonicalPolicy(source), executable_sha};
        require(Compact::Convert(source, binding, part, limits, std::cerr, false), "shard conversion failed");
      }
      else
      {
        auto* identity = file->GetDirectory(MetadataPath.c_str());
        require(identity, "compact identity missing");
        require(parameter<int>(*identity, "schema_version") == Compact::CollectionSchemaVersion &&
                    named(*identity, "schema_name") == Compact::CollectionSchemaName &&
                    named(*identity, "source_kind") == Compact::CollectionSourceKind,
                "not a converted shard collection");
        Compact::ExpectedSource binding{named(*identity, "source_file_sha256"), named(*identity, "source_configuration_fingerprint_fnv1a64"),
                                        named(*identity, "source_content_fingerprint_fnv1a64")};
        const auto source_manifest = named(*identity, "source_manifest_tsv");
        require(Compact::Sha256(source_manifest) == binding.file_sha256 && fingerprint(source_manifest) == binding.content_fingerprint_fnv1a64, "invalid source manifest binding");
        require(Compact::Read(file->GetDirectory("GarfieldFullTransport/Runtime"), geometry(*metadata), binding, part, limits, std::cerr, false), "invalid compact shard");
        require(Compact::ContentFingerprint(part) == input.content, "compact content differs from manifest");
      }
      require(fingerprint(named(*metadata, "canonical_configuration")) == part.source.configuration_fingerprint_fnv1a64, "physical fingerprint mismatch");
      if (expected.empty())
      {
        expected = common;
        data.geometry = part.geometry;
        data.source = part.source;
        data.source.file_sha256 = input_sha;
        data.source.content_fingerprint_fnv1a64 = fingerprint(manifest);
        data.source.converter_sha256 = executable_sha;
        data.first_root = data.end_root = part.first_root;
        data.max_depth = part.max_depth;
        data.max_nodes_per_root = part.max_nodes_per_root;
        data.max_queries_per_root = part.max_queries_per_root;
        data.guard_um = part.guard_um;
        data.root_lookup.assign(Map::RootCount(data.geometry), Compact::Invalid);
      }
      append_compact(data, part);
      for (unsigned int i = 0; i < counters.size(); ++i)
      {
        counters[i] += count(*metadata, GenerationCounters[i]);
      }
      unchanged(input);
    }
    require(data.nodes.size() == total_nodes && data.corners.size() == total_corners && data.values.size() <= total_values, "preflight/actual array counts differ");
    data.proof.sha256 = Compact::Sha256(Compact::ProofCanonical(data));
    require(Compact::Validate(data, limits, std::cerr, complete), "invalid collection or incomplete requested root coverage");
    const auto content = Compact::ContentFingerprint(data);
    const auto proof = data.proof;
    const auto nodes = data.nodes.size();
    const auto corners = data.corners.size();
    const auto values = data.values.size();
    const auto first = data.first_root;
    const auto end = data.end_root;
    const auto geometry_value = data.geometry;
    const auto binding = data.source;
    const auto temp = output + ".partial." + std::to_string(getpid());
    {
      TFile file(temp.c_str(), "CREATE");
      require(!file.IsZombie(), "cannot create collection output");
      auto* top = file.mkdir("GarfieldFullTransport");
      auto* identity = top->mkdir("Metadata");
      identity->cd();
      TParameter<int>("schema_version", Compact::CollectionSchemaVersion).Write();
      TNamed("schema_name", Compact::CollectionSchemaName).Write();
      TParameter<int>("source_schema_version", 4).Write();
      TNamed("compact_layout", Compact::Layout).Write();
      TNamed("source_kind", Compact::CollectionSourceKind).Write();
      TNamed("source_manifest_tsv", manifest.c_str()).Write();
      TNamed("source_file_sha256", binding.file_sha256.c_str()).Write();
      TNamed("source_configuration_fingerprint_fnv1a64", binding.configuration_fingerprint_fnv1a64.c_str()).Write();
      TNamed("source_content_fingerprint_fnv1a64", binding.content_fingerprint_fnv1a64.c_str()).Write();
      auto* original = top->mkdir("SourceMetadata");
      original->cd();
      for (const auto& [key, object] : saved)
      {
        (void) key;
        require(object->Write() > 0, "common metadata write failed");
      }
      TNamed("schema_name", Compact::CollectionSourceName).Write();
      TParameter<int>("generation_phi_chunk_index", 0).Write();
      TParameter<int>("generation_phi_chunk_count", 1).Write();
      TParameter<Long64_t>("generation_root_first", first).Write();
      TParameter<Long64_t>("generation_root_end", end).Write();
      TParameter<double>("maximum_drift_time_ns", data.maximum_drift_time_ns).Write();
      TNamed("reach_probe_trajectories", std::to_string(counters[0]).c_str()).Write();
      for (unsigned int i = 0; i < counters.size(); ++i)
      {
        TParameter<Long64_t>(GenerationCounters[i], counters[i]).Write();
      }
      require(Compact::Write(top->mkdir("Runtime"), data, limits, std::cerr), "compact array write failed");
      require(PHGarfieldFlushAndCloseRootOutput(file, std::cerr), "ROOT flush/close failed");
    }
    data = {};  // Never overlap owning arrays with independent ROOT readback.
    {
      auto file = open(temp);
      Compact::ExpectedSource expected_source{binding.file_sha256, binding.configuration_fingerprint_fnv1a64, binding.content_fingerprint_fnv1a64};
      require(Compact::Read(file->GetDirectory("GarfieldFullTransport/Runtime"), geometry_value, expected_source, data, limits, std::cerr, complete), "independent ROOT readback failed");
      require(Compact::ContentFingerprint(data) == content, "ROOT readback changed content");
    }
    require(Compact::FileSha256(get("--input-list")) == input_sha, "input manifest changed");
    for (const auto& input : inputs)
    {
      unchanged(input);
    }
    const auto output_sha = Compact::FileSha256(temp);
    require(::link(temp.c_str(), output.c_str()) == 0, "exclusive publication failed; partial preserved");
    require(::unlink(temp.c_str()) == 0, "partial unlink failed");
    std::ostringstream json;
    json << "{\"status\":\"PASS_EXACT_COMPACT_COLLECTION\",\"complete\":" << (complete ? "true" : "false")
         << ",\"root\":{\"path\":" << quote(std::filesystem::canonical(output)) << ",\"sha256\":" << quote(output_sha) << "}"
         << ",\"source_sha256\":" << quote(binding.file_sha256) << ",\"source_configuration\":" << quote(binding.configuration_fingerprint_fnv1a64)
         << ",\"source_content\":" << quote(binding.content_fingerprint_fnv1a64) << ",\"content_fingerprint\":" << quote(content)
         << ",\"first_root\":" << first << ",\"end_root\":" << end << ",\"source_points\":" << proof.source_points << ",\"nodes\":" << nodes
         << ",\"corners\":" << corners << ",\"values\":" << values << ",\"source_files\":" << inputs.size() << "}\n";
    const int fd = ::open(summary.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
    require(fd >= 0, "cannot publish exclusive summary");
    const auto text = json.str();
    require(::write(fd, text.data(), text.size()) == static_cast<ssize_t>(text.size()), "summary write failed");
    require(::fsync(fd) == 0 && ::close(fd) == 0, "summary close failed");
    std::cout << text;
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}

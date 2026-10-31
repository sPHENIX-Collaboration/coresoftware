#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TNamed.h>
#include <TParameter.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include "PHGarfieldCompactMap.h"
#include "PHGarfieldMapTools.h"
#include "PHGarfieldRootOutput.h"

namespace
{
  using Map = PHGarfieldCompactMap;
  using PHGarfieldMapTools::integer;
  using PHGarfieldMapTools::named;
  using PHGarfieldMapTools::parameter;
  using PHGarfieldMapTools::quote;
  using PHGarfieldMapTools::require;
  using PHGarfieldMapTools::write_exclusive;

  void copy_metadata(TDirectory* input, TDirectory* output)
  {
    TDirectory::TContext context(output);
    TIter next(input->GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(next()))
    {
      std::unique_ptr<TObject> object(key->ReadObj());
      require(object != nullptr, "unreadable source metadata");
      const std::string type = object->ClassName();
      require(type == "TNamed" || type.rfind("TParameter<", 0) == 0, "unexpected source metadata object " + type);
      require(object->Write(key->GetName()) > 0, "cannot copy source metadata");
    }
  }

}  // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--help")
  {
    std::cout << "Usage: " << argv[0] << " --input TRAINING.root --output COMPACT.root\n"
              << "  --proof-output CONVERSION.json --expect-source-sha256 SHA256\n"
              << "  --expect-source-configuration FNV1A64 --expect-source-content FNV1A64\n"
              << "\nConvert a complete adaptive schema4 map to compact schema5 arrays.\n"
              << "Preserves interpolation corner values and rejects existing outputs.\n"
              << "Optional unsigned resource ceilings:\n"
              << "  --max-source-points N --max-runtime-nodes N --max-runtime-values N\n"
              << "  --max-runtime-corners N --max-runtime-bytes N\n";
    return 0;
  }
  try
  {
    std::map<std::string, std::string> options;
    for (int i = 1; i < argc; i += 2)
    {
      require(i + 1 < argc && options.emplace(argv[i], argv[i + 1]).second, "options require unique key/value pairs");
    }
    const auto required = [&](const char* key) -> const std::string&
    {require(options.contains(key),std::string("missing ")+key);return options.at(key); };
    const std::vector<std::string> allowed{"--input", "--output", "--proof-output", "--expect-source-sha256", "--expect-source-configuration", "--expect-source-content",
                                           "--max-source-points", "--max-runtime-nodes", "--max-runtime-values", "--max-runtime-corners", "--max-runtime-bytes"};
    for (const auto& [key, value] : options)
    {
      (void) value;
      require(std::find(allowed.begin(), allowed.end(), key) != allowed.end(), "unsupported option " + key);
    }
    const auto input_path = std::filesystem::canonical(required("--input"));
    const auto output_path = std::filesystem::absolute(required("--output")).lexically_normal();
    const auto proof_path = std::filesystem::absolute(required("--proof-output")).lexically_normal();
    require(!std::filesystem::exists(output_path) && !std::filesystem::exists(proof_path), "preserve existing output/proof");
    require(output_path != proof_path, "ROOT and proof outputs must differ");
    Map::Limits limits;
    for (const auto& [key, target] : std::vector<std::pair<std::string, std::uint64_t*>>{
             {"--max-source-points", &limits.max_source_points}, {"--max-runtime-nodes", &limits.max_nodes}, {"--max-runtime-values", &limits.max_values}, {"--max-runtime-corners", &limits.max_corners}, {"--max-runtime-bytes", &limits.max_bytes}})
    {
      if (options.contains(key))
      {
        *target = integer(options.at(key));
      }
    }
    std::ostringstream errors;
    const bool valid_limits = Map::ValidateLimits(limits, errors);
    require(valid_limits, errors.str());
    Map::ExpectedSource expected{required("--expect-source-sha256"), required("--expect-source-configuration"), required("--expect-source-content")};
    require(Map::FileSha256(input_path) == expected.file_sha256, "input does not match externally pinned source SHA256");
    TFile input(input_path.c_str(), "READ");
    require(!input.IsZombie() && !input.TestBit(TFile::kRecovered), "training ROOT file unreadable/recovered");
    auto* metadata = input.GetDirectory("GarfieldFullTransport/Metadata");
    require(metadata, "source metadata missing");
    require(parameter<int>(*metadata, "schema_version") == 4 && named(*metadata, "schema_name") == PHGarfieldForwardMap::SchemaName, "converter requires complete schema4 training input");
    require(named(*metadata, "configuration_fingerprint_fnv1a64") == expected.configuration_fingerprint_fnv1a64, "source configuration differs from external expectation");
    const auto geometry = PHGarfieldMapTools::geometry(*metadata);
    PHGarfieldAdaptiveMap::Data source;
    require(PHGarfieldAdaptiveMap::Read(input.GetDirectory("GarfieldFullTransport/Adaptive"), geometry, source, errors, true), errors.str());
    require(PHGarfieldAdaptiveMap::ContentFingerprint(source) == expected.content_fingerprint_fnv1a64, "source content differs from external expectation");
    Map::SourceBinding binding{expected.file_sha256, expected.configuration_fingerprint_fnv1a64, expected.content_fingerprint_fnv1a64,
                               PHGarfieldAdaptiveMap::BuildCanonicalPolicy(source), Map::FileSha256("/proc/self/exe")};
    Map::Data data;
    require(Map::Convert(source, binding, data, limits, errors, true), errors.str());
    const auto temp_path = output_path.string() + ".temporary-" + std::to_string(getpid());
    {
      TFile output(temp_path.c_str(), "CREATE");
      require(!output.IsZombie(), "cannot create unique temporary ROOT output");
      auto* top = output.mkdir("GarfieldFullTransport");
      require(top, "cannot create transport directory");
      auto* identity = top->mkdir("Metadata");
      require(identity, "cannot create compact identity metadata");
      {
        TDirectory::TContext context(identity);
        require(TParameter<int>("schema_version", Map::SchemaVersion).Write() > 0 && TNamed("schema_name", Map::SchemaName).Write() > 0 &&
                    TParameter<int>("source_schema_version", 4).Write() > 0 && TNamed("compact_layout", Map::Layout).Write() > 0 &&
                    TNamed("source_file_sha256", binding.file_sha256.c_str()).Write() > 0 &&
                    TNamed("source_configuration_fingerprint_fnv1a64", binding.configuration_fingerprint_fnv1a64.c_str()).Write() > 0 &&
                    TNamed("source_content_fingerprint_fnv1a64", binding.content_fingerprint_fnv1a64.c_str()).Write() > 0,
                "cannot write compact identity metadata");
      }
      auto* copied = top->mkdir("SourceMetadata");
      require(copied, "cannot create source metadata directory");
      copy_metadata(metadata, copied);
      require(Map::Write(top->mkdir("Runtime"), data, limits, errors), errors.str());
      require(PHGarfieldFlushAndCloseRootOutput(output, errors), errors.str());
    }
    const auto content = Map::ContentFingerprint(data);
    const auto proof = data.proof;
    const auto nodes = data.nodes.size();
    const auto corners = data.corners.size();
    const auto values = data.values.size();
    data = {};  // Release converted vectors before independent ROOT readback.
    {
      TFile check(temp_path.c_str(), "READ");
      require(!check.IsZombie() && !check.TestBit(TFile::kRecovered), "compact ROOT output reopened recovered/unreadable");
      require(Map::Read(check.GetDirectory("GarfieldFullTransport/Runtime"), geometry, expected, data, limits, errors, true), errors.str());
      require(Map::ContentFingerprint(data) == content && Map::VerifyConversion(source, data, errors), errors.str());
    }
    input.Close();
    require(Map::FileSha256(input_path) == expected.file_sha256, "source changed during conversion");
    const auto output_sha = Map::FileSha256(temp_path);
    require(::link(temp_path.c_str(), output_path.c_str()) == 0, "atomic ROOT publication failed; preserved temporary file");
    require(::unlink(temp_path.c_str()) == 0, "published ROOT but temporary unlink failed");
    std::ostringstream json;
    json << "{\n  \"schema\":\"phgarfield_compact_conversion_manifest_v1\",\n  \"status\":\"PASS_CONVERSION_ONLY\",\n"
         << "  \"source\":{\"path\":" << quote(input_path) << ",\"sha256\":" << quote(expected.file_sha256) << "},\n"
         << "  \"compact\":{\"path\":" << quote(output_path) << ",\"sha256\":" << quote(output_sha) << "},\n"
         << "  \"source_configuration_fingerprint_fnv1a64\":" << quote(expected.configuration_fingerprint_fnv1a64) << ",\n"
         << "  \"source_content_fingerprint_fnv1a64\":" << quote(expected.content_fingerprint_fnv1a64) << ",\n"
         << "  \"compact_content_fingerprint_fnv1a64\":" << quote(content) << ",\n  \"conversion_proof_sha256\":" << quote(proof.sha256)
         << ",\n  \"converter_sha256\":" << quote(binding.converter_sha256) << ",\n  \"nodes\":" << nodes << ",\n  \"certified_corner_records\":" << corners << ",\n  \"unique_corner_values\":" << values << ",\n"
         << "  \"runtime_limits\":{\"max_source_points\":" << limits.max_source_points << ",\"max_nodes\":" << limits.max_nodes
         << ",\"max_corners\":" << limits.max_corners << ",\"max_values\":" << limits.max_values << ",\"max_roots\":" << limits.max_roots << ",\"max_bytes\":" << limits.max_bytes << "},\n"
         << "  \"checks\":{\"all_source_bounds_routes_corner_bits_exact\":true,\"independent_ROOT_readback_unrecovered\":true,\"source_immutable\":true},\n"
         << "  \"scope\":\"Conversion and structure only; source accuracy/coverage and detector-production approval remain external prerequisites.\"\n}\n";
    const auto temporary_proof = proof_path.string() + ".temporary-" + std::to_string(getpid());
    write_exclusive(temporary_proof, json.str());
    require(::link(temporary_proof.c_str(), proof_path.c_str()) == 0, "atomic proof publication failed; ROOT/proof evidence preserved");
    require(::unlink(temporary_proof.c_str()) == 0, "proof temporary unlink failed");
    std::cout << json.str();
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "Compact conversion failed: " << error.what() << '\n';
    return 1;
  }
}

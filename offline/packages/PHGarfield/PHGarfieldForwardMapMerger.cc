#include "PHGarfieldForwardMap.h"
#include "PHGarfieldMapTools.h"
#include "PHGarfieldRootOutput.h"

#include <TAxis.h>
#include <TH3.h>
#include <TNamed.h>
#include <TParameter.h>

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  using PHGarfieldMapTools::common_metadata;
  using PHGarfieldMapTools::fingerprint;
  using PHGarfieldMapTools::integer;
  using PHGarfieldMapTools::MetadataPath;
  using PHGarfieldMapTools::named;
  using PHGarfieldMapTools::Objects;
  using PHGarfieldMapTools::open;
  using PHGarfieldMapTools::parameter;
  using PHGarfieldMapTools::require;

  struct Shard
  {
    std::string path;
    std::unique_ptr<TFile> file;
    int index{0};
  };

  bool same_axis(const TAxis& first, const TAxis& second)
  {
    if (first.GetNbins() != second.GetNbins())
    {
      return false;
    }
    for (int bin = 1; bin <= first.GetNbins() + 1; ++bin)
    {
      if (first.GetBinLowEdge(bin) != second.GetBinLowEdge(bin))
      {
        return false;
      }
    }
    return true;
  }

  void merge(const std::string& destination, const std::vector<std::string>& paths)
  {
    require(!std::filesystem::exists(destination), "output already exists");
    std::vector<Shard> shards;
    std::map<std::string, std::string> expected;
    Objects saved;
    int count_expected = 0;
    int phi_nodes = 0;
    int version = 0;
    double maximum_time = 0.;
    std::uint64_t probes = 0;
    std::map<int, std::string> identities;
    for (const auto& path : paths)
    {
      auto file = open(path);
      auto* metadata = file->GetDirectory(MetadataPath.c_str());
      require(metadata, "missing forward-map metadata");
      const int schema = parameter<int>(*metadata, "schema_version");
      require(schema >= 1 && schema <= 3 && named(*metadata, "schema_name") == "PHG4TpcGarfieldFullTransportLUTShard", "requires regular-grid shards (schemas 1-3)");
      const int index = parameter<int>(*metadata, "generation_phi_chunk_index");
      const int chunks = parameter<int>(*metadata, "generation_phi_chunk_count");
      const int phi = parameter<int>(*metadata, "phi_nodes");
      require(chunks > 0 && chunks <= phi && phi <= std::numeric_limits<int>::max() - 2 && index >= 0 && index < chunks && identities.emplace(index, path).second, "invalid or duplicate shard identity");
      const auto common = common_metadata(*metadata, expected.empty() ? &saved : nullptr);
      if (expected.empty())
      {
        expected = common;
        count_expected = chunks;
        phi_nodes = phi;
        version = schema;
      }
      require(expected == common && chunks == count_expected && schema == version, "incompatible shard geometry, physics or schema");
      const auto canonical = named(*metadata, "canonical_configuration");
      require(fingerprint(canonical) == named(*metadata, "configuration_fingerprint_fnv1a64"), "configuration fingerprint mismatch");
      const double duration = parameter<double>(*metadata, "maximum_drift_time_ns");
      require(std::isfinite(duration) && duration >= 0., "invalid maximum drift time");
      maximum_time = std::max(maximum_time, duration);
      const auto amount = integer(named(*metadata, "reach_probe_trajectories"));
      require(amount <= std::numeric_limits<std::uint64_t>::max() - probes, "probe count overflow");
      probes += amount;
      shards.push_back({path, std::move(file), index});
    }
    require(shards.size() == static_cast<std::size_t>(count_expected), "missing shards");
    std::sort(shards.begin(), shards.end(), [](const auto& a, const auto& b)
              { return a.index < b.index; });
    const auto temporary = destination + ".partial." + std::to_string(getpid());
    {
      TFile output(temporary.c_str(), "CREATE");
      require(!output.IsZombie(), "cannot create temporary ROOT output");
      auto* transport = output.mkdir("GarfieldFullTransport");
      require(transport, "cannot create transport directory");
      std::vector<std::string> stems{"hDeltaR", "hDeltaRPhi", "hDriftTime", "hReachesReadoutNode", "hReachCellClass"};
      if (version >= 3)
      {
        stems.insert(stems.end(), {"hInterpolationTransverseError", "hInterpolationTimeEquivalentError"});
      }
      for (const std::string side : {"negz", "posz"})
      {
        for (const auto& stem : stems)
        {
          const auto name = stem + "_" + side;
          const auto key = "GarfieldFullTransport/" + name;
          auto* first = dynamic_cast<TH3*>(shards.front().file->Get(key.c_str()));
          require(first, "missing histogram " + key);
          std::unique_ptr<TH3> histogram(dynamic_cast<TH3*>(first->Clone(name.c_str())));
          require(histogram != nullptr, "cannot clone histogram");
          histogram->SetDirectory(nullptr);
          histogram->Reset();
          const bool cells = stem == "hReachCellClass" || stem.starts_with("hInterpolation");
          require(first->GetNbinsX() == phi_nodes + (cells ? 0 : 2), "invalid phi support dimensions");
          for (const auto& shard : shards)
          {
            auto* part = dynamic_cast<TH3*>(shard.file->Get(key.c_str()));
            require(part && std::string(part->ClassName()) == first->ClassName() && same_axis(*first->GetXaxis(), *part->GetXaxis()) && same_axis(*first->GetYaxis(), *part->GetYaxis()) && same_axis(*first->GetZaxis(), *part->GetZaxis()), "incompatible histogram " + key);
            const int begin = static_cast<int>(static_cast<std::uint64_t>(phi_nodes) * static_cast<std::uint64_t>(shard.index) / static_cast<std::uint64_t>(count_expected)) + (cells ? 1 : 2);
            const int end = static_cast<int>(static_cast<std::uint64_t>(phi_nodes) * static_cast<std::uint64_t>(shard.index + 1) / static_cast<std::uint64_t>(count_expected)) + (cells ? 1 : 2);
            for (int x = begin; x < end; ++x)
            {
              for (int y = 1; y <= histogram->GetNbinsY(); ++y)
              {
                for (int z = 1; z <= histogram->GetNbinsZ(); ++z)
                {
                  histogram->SetBinContent(x, y, z, part->GetBinContent(x, y, z));
                }
              }
            }
          }
          if (!cells)
          {
            const int nx = histogram->GetNbinsX();
            for (int y = 1; y <= histogram->GetNbinsY(); ++y)
            {
              for (int z = 1; z <= histogram->GetNbinsZ(); ++z)
              {
                histogram->SetBinContent(1, y, z, histogram->GetBinContent(nx - 1, y, z));
                histogram->SetBinContent(nx, y, z, histogram->GetBinContent(2, y, z));
              }
            }
          }
          transport->cd();
          require(histogram->Write() > 0, "histogram write failed");
        }
      }
      auto* metadata = transport->mkdir("Metadata");
      require(metadata, "cannot create metadata directory");
      metadata->cd();
      for (const auto& [key, object] : saved)
      {
        (void) key;
        require(object->Write() > 0, "metadata write failed");
      }
      require(TNamed("schema_name", PHGarfieldForwardMap::SchemaName).Write() > 0 && TParameter<int>("generation_phi_chunk_index", 0).Write() > 0 && TParameter<int>("generation_phi_chunk_count", 1).Write() > 0 && TParameter<double>("maximum_drift_time_ns", maximum_time).Write() > 0 && TNamed("reach_probe_trajectories", std::to_string(probes).c_str()).Write() > 0, "aggregate metadata write failed");
      require(PHGarfieldFlushAndCloseRootOutput(output, std::cerr), "ROOT write/close failed");
    }
    auto check = open(temporary);
    require(PHGarfieldForwardMap::ValidateRegular(*check, std::cerr), "invalid merged regular map");
    check->Close();
    require(::link(temporary.c_str(), destination.c_str()) == 0, "exclusive publication failed");
    require(::unlink(temporary.c_str()) == 0, "temporary output unlink failed");
    std::cout << "Merged " << shards.size() << " regular-grid shards into " << destination << '\n';
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--help")
  {
    std::cout << "Usage: PHGarfieldForwardMapMerger OUTPUT.root SHARD.root [SHARD.root ...]\n";
    return 0;
  }
  try
  {
    PHGarfieldMapTools::require(argc >= 3, "output and at least one source shard required");
    merge(argv[1], std::vector<std::string>(argv + 2, argv + argc));
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}

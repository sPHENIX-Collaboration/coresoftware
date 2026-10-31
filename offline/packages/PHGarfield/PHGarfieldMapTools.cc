#include "PHGarfieldMapTools.h"
#include <cerrno>

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

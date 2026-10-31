#ifndef PHGARFIELDMAPTOOLS_H
#define PHGARFIELDMAPTOOLS_H

// Internal helpers shared by the standalone forward-map conversion tools.
#include "PHGarfieldAdaptiveMap.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TObject.h>
#include <TParameter.h>

#include <array>
#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>

namespace PHGarfieldMapTools
{
  using Map = PHGarfieldAdaptiveMap;
  using Objects = std::map<std::string, std::unique_ptr<TObject>>;
  extern const std::string MetadataPath;
  extern const std::string AdaptivePath;
  extern const std::array<const char*, 5> GenerationCounters;

  struct Capacity
  {
    bool enabled{false};
    std::uint64_t points{0}, cells{0}, bytes{0};
    std::string proof_output;
  };

  void require(bool value, const std::string& message);
  template <class T>
  T parameter(TDirectory& directory, const char* name)
  {
    const auto* object = dynamic_cast<TParameter<T>*>(directory.Get(name));
    require(object != nullptr, std::string("missing/incorrect metadata ") + name);
    return object->GetVal();
  }

  std::uint64_t count(TDirectory& directory, const char* name);
  std::string named(TDirectory& directory, const char* name);
  std::string fingerprint(const std::string& text);
  std::uint64_t integer(const std::string& text);
  void apply_capacity(Map::Policy& policy, const Capacity& capacity);
  std::string capacity_prefix(const std::string& canonical);
  std::array<std::uint64_t, 3> canonical_capacity(const std::string& canonical);
  std::string primitive_fingerprint(const Map::Data& data);
  bool same_values(const std::array<double, 3>& first, const std::array<double, 3>& second);
  std::string quote(const std::string& value);
  void write_exclusive(const std::string& path, const std::string& text);
  void check_capacity_envelope(TDirectory& metadata, const Map::Data& data);
  Map::Geometry geometry(TDirectory& directory);
  std::unique_ptr<TFile> open(const std::string& path);
  void check_summary(TDirectory& metadata, const Map::Data& data);
  std::array<std::uint64_t, 5> generation_counters(TDirectory& metadata, const Map::Data& data);
  void check_file(const std::string& path, bool complete);
  std::string signature(const TObject& object);
  std::map<std::string, std::string> common_metadata(TDirectory& directory, Objects* saved);
}  // namespace PHGarfieldMapTools

#endif

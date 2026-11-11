#ifndef PHGARFIELDMAPTOOLS_H
#define PHGARFIELDMAPTOOLS_H

// Internal helpers shared by the standalone forward-map conversion tools.
#include "PHGarfieldCompactMap.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TObject.h>
#include <TParameter.h>

#include <array>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>

namespace PHGarfieldMapTools
{
  using Map = PHGarfieldAdaptiveMap;
  using Objects = std::map<std::string, std::unique_ptr<TObject>>;
  using Options = std::map<std::string, std::string>;
  extern const std::string MetadataPath;
  extern const std::string AdaptivePath;
  extern const std::array<const char*, 5> GenerationCounters;

  void require(bool value, const std::string& message);
  Options options(int argc, char** argv, std::initializer_list<std::string> allowed);
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
  std::string quote(const std::string& value);
  void write_exclusive(const std::string& path, const std::string& text);
  void publish_exclusive(const std::string& temporary, const std::string& output);
  void write_compact_output(const std::string& path, const PHGarfieldCompactMap::Data& data,
                            const PHGarfieldCompactMap::Limits& limits, std::ostream& errors,
                            const std::function<void(TDirectory&)>& write_source_metadata,
                            const std::string* collection_manifest = nullptr);
  // Release the owning arrays before independently reopening and checking them.
  void readback_compact_output(const std::string& path, PHGarfieldCompactMap::Data& data,
                               const PHGarfieldCompactMap::Limits& limits, std::ostream& errors,
                               bool complete);
  Map::Geometry geometry(TDirectory& directory);
  std::unique_ptr<TFile> open(const std::string& path);
  void check_summary(TDirectory& metadata, const Map::Data& data);
  std::string signature(const TObject& object);
  std::map<std::string, std::string> common_metadata(TDirectory& directory, Objects* saved);
}  // namespace PHGarfieldMapTools

#endif

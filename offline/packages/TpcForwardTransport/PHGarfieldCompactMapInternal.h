#ifndef PHGARFIELD_COMPACTMAPINTERNAL_H
#define PHGARFIELD_COMPACTMAPINTERNAL_H

#include "PHGarfieldCompactMap.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <ostream>
#include <string>

// Shared bit-layout and exact source checks for runtime validation and tools.
namespace PHGarfieldCompactMapInternal
{
  using Map = PHGarfieldCompactMap;
  inline bool require(bool value, std::ostream& error, const char* text)
  {
    if (!value)
    {
      error << "Compact LUT: " << text << '\n';
    }
    return value;
  }
  inline bool hex(const std::string& value, std::size_t length)
  {
    return value.size() == length && std::all_of(value.begin(), value.end(), [](char c)
                                                 { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); });
  }
  inline bool bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
  inline bool bounds_equal(const Map::Bounds& a, const Map::Bounds& b)
  {
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      if (!bits(a.lo[axis], b.lo[axis]) || !bits(a.hi[axis], b.hi[axis]))
      {
        return false;
      }
    }
    return true;
  }
  inline bool geometry_equal(const Map::Geometry& a, const Map::Geometry& b)
  {
    return a.phi_nodes == b.phi_nodes && a.radial_nodes == b.radial_nodes && a.longitudinal_nodes == b.longitudinal_nodes &&
           bits(a.inner_radius_cm, b.inner_radius_cm) && bits(a.outer_radius_cm, b.outer_radius_cm) &&
           bits(a.central_membrane_half_width_cm, b.central_membrane_half_width_cm) && bits(a.readout_abs_z_cm, b.readout_abs_z_cm) &&
           bits(a.nominal_velocity_cm_ns, b.nominal_velocity_cm_ns);
  }
  inline bool binding_valid(const Map::SourceBinding& b)
  {
    return hex(b.file_sha256, 64) && hex(b.configuration_fingerprint_fnv1a64, 16) && hex(b.content_fingerprint_fnv1a64, 16) &&
           hex(b.converter_sha256, 64) && !b.canonical_policy.empty() && b.canonical_policy.size() <= 65536;
  }
  inline unsigned int kind(const Map::Node& n) { return n.tag & 3U; }
  inline int leaf_class(const Map::Node& n) { return static_cast<int>((n.tag >> 2U) & 7U); }
  inline int reason(const Map::Node& n) { return static_cast<int>((n.tag >> 5U) & 7U); }
  inline bool class_reason(int c, int r)
  {
    if (c == 0)
    {
      return r == 5;
    }
    if (c == 1)
    {
      return r == 0;
    }
    if (c == 2)
    {
      return r == 1 || r == 3 || r == 4;
    }
    if (c == 3)
    {
      return r == 7;
    }
    if (c == 4)
    {
      return r == 2 || r == 3 || r == 4 || r == 6;
    }
    return false;
  }
}  // namespace PHGarfieldCompactMapInternal

#endif

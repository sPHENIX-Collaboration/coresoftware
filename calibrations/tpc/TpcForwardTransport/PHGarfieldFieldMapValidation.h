#ifndef PHGARFIELDFIELDMAPVALIDATION_H
#define PHGARFIELDFIELDMAPVALIDATION_H

#include <array>
#include <iosfwd>
#include <string>

// Basic inspection of the Cartesian ROOT map consumed by PHField3DCartesian.
// This helper inspects a filename; field ownership and CDB resolution stay with
// PHGarfieldForwardTransport.
namespace PHGarfieldFieldMapValidation
{
  struct CartesianDomain
  {
    std::array<double, 3> minimum_cm{};
    std::array<double, 3> maximum_cm{};
  };

  // Bounds are the inclusive coordinate extrema in the map frame, in cm.
  // Requires readable finite entries and a nondegenerate range on every axis.
  // The input's grid spacing, completeness and uniqueness are trusted.
  // On failure both arrays contain NaN, including when a valid result is reused.
  bool InspectCartesianMap(const std::string& filename,
                           CartesianDomain& domain,
                           std::ostream& diagnostics);
}  // namespace PHGarfieldFieldMapValidation

#endif

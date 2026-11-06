#ifndef PHGARFIELDFIELDMAPVALIDATION_H
#define PHGARFIELDFIELDMAPVALIDATION_H

#include <array>
#include <iosfwd>
#include <string>

// Internal validation of the Cartesian ROOT grid consumed by PHField3DCartesian.
// This helper inspects a filename; field ownership and CDB resolution stay with
// PHGarfield.
namespace PHGarfieldFieldMapValidation
{
  struct CartesianDomain
  {
    std::array<double, 3> minimum_cm{};
    std::array<double, 3> maximum_cm{};
  };

  // Bounds are the inclusive coordinate extrema in the map frame, in cm.
  // On failure both arrays contain NaN, including when a valid result is reused.
  bool InspectCartesianMap(const std::string& filename,
                           CartesianDomain& domain,
                           std::ostream& diagnostics);
}  // namespace PHGarfieldFieldMapValidation

#endif

#include "PHGarfieldFieldMapValidation.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <vector>

bool PHGarfieldFieldMapValidation::InspectCartesianMap(
    const std::string& filename,
    CartesianDomain& domain,
    std::ostream& diagnostics)
{
  const double invalid = std::numeric_limits<double>::quiet_NaN();
  domain.minimum_cm.fill(invalid);
  domain.maximum_cm.fill(invalid);
  if (filename.empty())
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap - Magnetic-field map filename is empty"
                << std::endl;
    return false;
  }

  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap - Could not inspect magnetic-field map: "
                << filename << std::endl;
    return false;
  }

  auto* field_map = dynamic_cast<TTree*>(input->Get("fieldmap"));
  if (!field_map)
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                << " Magnetic-field map has no fieldmap tree: "
                << filename << std::endl;
    return false;
  }

  const std::array<const char*, 6> branch_names =
      {{"x", "y", "z", "bx", "by", "bz"}};
  for (const char* branch_name : branch_names)
  {
    if (!field_map->GetLeaf(branch_name))
    {
      diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap - Magnetic-field map is missing branch "
                  << branch_name << ": " << filename
                  << std::endl;
      return false;
    }
  }

  Float_t x_cm = 0.0F;
  Float_t y_cm = 0.0F;
  Float_t z_cm = 0.0F;
  Float_t bx_t = 0.0F;
  Float_t by_t = 0.0F;
  Float_t bz_t = 0.0F;
  if (field_map->SetBranchAddress("x", &x_cm) < 0 ||
      field_map->SetBranchAddress("y", &y_cm) < 0 ||
      field_map->SetBranchAddress("z", &z_cm) < 0 ||
      field_map->SetBranchAddress("bx", &bx_t) < 0 ||
      field_map->SetBranchAddress("by", &by_t) < 0 ||
      field_map->SetBranchAddress("bz", &bz_t) < 0)
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                << " Failed to bind magnetic-field map branches: "
                << filename << std::endl;
    return false;
  }

  const Long64_t entries = field_map->GetEntries();
  constexpr Long64_t maximum_supported_entries = 100000000;
  if (entries <= 0 || entries > maximum_supported_entries)
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                << " Unsupported magnetic-field grid size " << entries
                << ": " << filename << std::endl;
    return false;
  }

  std::array<std::set<float>, 3> axis_values;
  std::vector<std::array<float, 3>> coordinates;
  coordinates.reserve(static_cast<std::size_t>(entries));
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    if (field_map->GetEntry(entry) <= 0 ||
        !std::isfinite(x_cm) ||
        !std::isfinite(y_cm) ||
        !std::isfinite(z_cm) ||
        !std::isfinite(bx_t) ||
        !std::isfinite(by_t) ||
        !std::isfinite(bz_t))
    {
      diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                  << " Invalid magnetic-field grid entry " << entry
                  << ": " << filename
                  << std::endl;
      return false;
    }
    axis_values[0].insert(x_cm);
    axis_values[1].insert(y_cm);
    axis_values[2].insert(z_cm);
    coordinates.push_back({{x_cm, y_cm, z_cm}});
  }

  const auto is_uniform_axis =
      [](const std::set<float>& values)
  {
    if (values.size() < 2)
    {
      return false;
    }
    const double expected_spacing =
        (static_cast<double>(*values.rbegin()) -
         static_cast<double>(*values.begin())) /
        static_cast<double>(values.size() - 1);
    auto previous = values.begin();
    for (auto current = std::next(previous);
         current != values.end();
         ++current)
    {
      const double spacing =
          static_cast<double>(*current) -
          static_cast<double>(*previous);
      const double scale =
          std::max(1.0, std::abs(expected_spacing));
      if (!std::isfinite(spacing) ||
          std::abs(spacing - expected_spacing) > 1.0e-6 * scale)
      {
        return false;
      }
      previous = current;
    }
    return true;
  };
  if (!std::all_of(
          axis_values.begin(),
          axis_values.end(),
          is_uniform_axis))
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                << " Magnetic-field coordinates do not form uniform axes: "
                << filename << std::endl;
    return false;
  }

  std::size_t expected_entries = 1;
  for (const auto& values : axis_values)
  {
    if (expected_entries >
        std::numeric_limits<std::size_t>::max() / values.size())
    {
      diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                  << " Magnetic-field grid dimensions overflow: "
                  << filename << std::endl;
      return false;
    }
    expected_entries *= values.size();
  }
  if (expected_entries != coordinates.size())
  {
    diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                << " Magnetic-field map is not a complete Cartesian lattice: "
                << filename << std::endl;
    return false;
  }

  std::array<std::map<float, std::size_t>, 3> axis_index;
  for (unsigned int axis = 0; axis < axis_values.size(); ++axis)
  {
    std::size_t index = 0;
    for (const float value : axis_values[axis])
    {
      axis_index[axis].emplace(value, index++);
    }
  }

  std::vector<unsigned char> occupied(expected_entries, 0);
  const std::size_t ny = axis_values[1].size();
  const std::size_t nz = axis_values[2].size();
  for (const auto& coordinate : coordinates)
  {
    const std::size_t linear_index =
        (axis_index[0].at(coordinate[0]) * ny +
         axis_index[1].at(coordinate[1])) *
            nz +
        axis_index[2].at(coordinate[2]);
    if (occupied[linear_index] != 0)
    {
      diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                  << " Magnetic-field map contains a duplicate grid point: "
                  << filename << std::endl;
      return false;
    }
    occupied[linear_index] = 1;
  }

  // Publish the domain only after the full lattice has been validated.
  for (unsigned int axis = 0; axis < axis_values.size(); ++axis)
  {
    domain.minimum_cm[axis] =
        static_cast<double>(*axis_values[axis].begin());
    domain.maximum_cm[axis] =
        static_cast<double>(*axis_values[axis].rbegin());
  }
  return true;
}

#include "PHGarfieldFieldMapValidation.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <ostream>

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

  CartesianDomain inspected;
  inspected.minimum_cm.fill(std::numeric_limits<double>::max());
  inspected.maximum_cm.fill(std::numeric_limits<double>::lowest());
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
    const std::array<double, 3> coordinate = {{x_cm, y_cm, z_cm}};
    for (std::size_t axis = 0; axis < coordinate.size(); ++axis)
    {
      inspected.minimum_cm[axis] = std::min(inspected.minimum_cm[axis], coordinate[axis]);
      inspected.maximum_cm[axis] = std::max(inspected.maximum_cm[axis], coordinate[axis]);
    }
  }

  for (std::size_t axis = 0; axis < inspected.minimum_cm.size(); ++axis)
  {
    if (inspected.minimum_cm[axis] >= inspected.maximum_cm[axis])
    {
      diagnostics << "PHGarfieldFieldMapValidation::InspectCartesianMap -"
                  << " Magnetic-field map has a degenerate coordinate range on axis "
                  << axis << ": " << filename << std::endl;
      return false;
    }
  }

  // Trust the supplied map's grid structure; coverage checks use its extrema.
  // Do not audit grid spacing, missing Cartesian points or duplicate entries.
  domain = inspected;
  return true;
}

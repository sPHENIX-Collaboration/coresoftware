#include "PHG4TpcGeom.h"

namespace
{
  const std::array<std::vector<double>, 2> dummy_array;
}

const std::array<std::vector<double>, 2> &PHG4TpcGeom::get_sector_min_phi()
  {
    PHOOL_VIRTUAL_WARN("get_sector_min_phi()");
    return dummy_array;
  }

const std::array<std::vector<double>, 2> &PHG4TpcGeom::get_sector_max_phi()
  {
    PHOOL_VIRTUAL_WARN("get_sector_max_phi()");
    return dummy_array;
  }

void PHG4TpcGeom::identify(std::ostream& os) const
{
  os << "virtual PHG4TpcGeom"
     << std::endl;
  return;
}

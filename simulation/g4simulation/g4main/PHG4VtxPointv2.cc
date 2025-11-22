#include "PHG4VtxPointv2.h"

void PHG4VtxPointv2::identify(std::ostream& os) const
{
  os << "vtx position: x: " << get_x()
     << ", y: " << get_y()
     << ", z: " << get_z()
     << ", t: " << get_t()
     << ", prodProcess: " << getProdProcessAsString().data()
     << std::endl;
}

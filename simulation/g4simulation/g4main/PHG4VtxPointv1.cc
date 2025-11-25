#include "PHG4VtxPointv1.h"

void PHG4VtxPointv1::identify(std::ostream& os) const
{
  os << "vtx position: x: " << get_x()
     << ", y: " << get_y()
     << ", z: " << get_z()
     << ", t: " << get_t()
     << std::endl;
}

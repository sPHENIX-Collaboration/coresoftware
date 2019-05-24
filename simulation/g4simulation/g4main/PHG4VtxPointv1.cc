#include "PHG4VtxPointv1.h"

#include <ostream>  // for operator<<, basic_ostream, basic_ostream<>::__ost...

using namespace std;

void
PHG4VtxPointv1::identify(ostream& os) const
{
  os << "vtx position: x: " << get_x()
       << ", y: " << get_y()
       << ", z: " << get_z()
       << ", t: " << get_t()
     << endl;
}

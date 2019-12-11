#include "PHG4VtxPointv1.h"

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

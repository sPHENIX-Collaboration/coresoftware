#include "PHG4CylinderCell_MVTX.h"

using namespace std;

PHG4CylinderCell_MVTX::PHG4CylinderCell_MVTX()
: PHG4CylinderCellv2()
{}

void
PHG4CylinderCell_MVTX::identify(std::ostream& os) const
{
  os << "PHG4CylinderCell_MVTX: ";
  os << " layer, stave number, half stave number, module, chip, pixel =  "
     << layer << ","; 
  os << get_stave_index() << ",";
  os << get_half_stave_index() << ",";
  os << get_module_index() << ", ";
  os << get_chip_index() << ", ";
  os << get_pixel_index();
  os << endl;
}

#include "PHG4CylinderCell_MAPS.h"

using namespace std;

PHG4CylinderCell_MAPS::PHG4CylinderCell_MAPS()
: PHG4CylinderCellv2()
{}

void
PHG4CylinderCell_MAPS::identify(std::ostream& os) const
{
  os << "PHG4CylinderCell_MAPS: ";
  os << " layer, stave number, half stave number, module, chip, pixel =  "
     << layer << ","; 
  os << get_stave_index() << ",";
  os << get_half_stave_index() << ",";
  os << get_module_index() << ", ";
  os << get_chip_index() << ", ";
  os << get_pixel_index();
  os << endl;
}

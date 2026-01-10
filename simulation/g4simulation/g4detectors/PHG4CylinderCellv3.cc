#include "PHG4CylinderCellv3.h"

void PHG4CylinderCellv3::identify(std::ostream& os) const
{
  os << "PHG4CylinderCellv3: #" << cellid << " ";
  os << "(layer,e,j_index,k_index,l_index) = (";
  os << layer << ",";
  os << get_edep() << ",";
  os << get_j_index() << ",";
  os << get_k_index() << ",";
  os << get_l_index() << ",";
  os << ")";
  os << std::endl;
}

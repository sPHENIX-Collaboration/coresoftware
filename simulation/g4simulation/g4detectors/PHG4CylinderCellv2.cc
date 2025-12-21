#include "PHG4CylinderCellv2.h"

PHG4CylinderCellv2::PHG4CylinderCellv2()
  : PHG4CylinderCellv1()
{
}

void PHG4CylinderCellv2::identify(std::ostream& os) const
{
  os << "PHG4CylinderCellv2: #" << cellid << " ";
  os << "(layer,binz,binphi,e,sensor_index,phi_index,z_index) = (";
  os << layer << ",";
  os << binz << ",";
  os << binphi << ",";
  os << get_edep() << ",";
  os << get_sensor_index() << ",";
  os << get_ladder_phi_index() << ",";
  os << get_ladder_z_index();
  os << ")";
  os << std::endl;
}

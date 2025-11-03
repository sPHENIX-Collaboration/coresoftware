#include "TruthVertex_v1.h"

void TruthVertex_v1::identify(std::ostream& os) const
{  
  
  os << "---TruthVertex_v1--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << std::endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << std::endl;
  os << " x =  " << get_x() << " +/- " << get_x_err() << std::endl;
  os << " y =  " << get_y() << " +/- " << get_y_err() << std::endl;
  os << "-----------------------------------------------" << std::endl;

}

int TruthVertex_v1::isValid() const
{
  return std::isfinite(_z) && std::isfinite(_t);
}

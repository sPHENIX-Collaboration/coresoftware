#include "MbdVertexv2.h"

#include <cmath>

void MbdVertexv2::identify(std::ostream& os) const
{
  os << "---MbdVertexv2--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << std::endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int MbdVertexv2::isValid() const
{
  if (_id == std::numeric_limits<unsigned int>::max())
  {
    return 0;
  }
  if (std::isnan(_t))
  {
    return 0;
  }
  if (std::isnan(_t_err))
  {
    return 0;
  }
  if (std::isnan(_z))
  {
    return 0;
  }
  if (std::isnan(_z_err))
  {
    return 0;
  }

  return 1;
}

float MbdVertexv2::get_position(unsigned int coor) const
{
  if (coor == 0)
  {
    return get_x();
  }
  else if (coor == 1)
  {
    return get_y();
  }
  else if (coor == 2)
  {
    return get_z();
  }
  else
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
}

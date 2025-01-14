#include "MbdVertexv1.h"

#include <cmath>

void MbdVertexv1::identify(std::ostream& os) const
{
  os << "---MbdVertexv1--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << std::endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int MbdVertexv1::isValid() const
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

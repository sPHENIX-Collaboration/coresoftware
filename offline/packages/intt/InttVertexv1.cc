#include "InttVertexv1.h"

#include <cmath>

void InttVertexv1::identify(std::ostream& os) const
{
  os << "---InttVertexv1--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int InttVertexv1::isValid() const
{
  if (_id == std::numeric_limits<unsigned int>::max())
  {
    return 0;
  }
  if (isnan(_z))
  {
    return 0;
  }
  if (isnan(_z_err))
  {
    return 0;
  }

  return 1;
}

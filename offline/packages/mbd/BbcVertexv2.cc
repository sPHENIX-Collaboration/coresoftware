#include "BbcVertexv2.h"

#include <cmath>

BbcVertexv2::BbcVertexv2()
{
  for (int i = 0; i < 2; i++)
  {
    _bbc_ns_npmt[i] = 0;
    _bbc_ns_q[i] = std::numeric_limits<float>::quiet_NaN();
    _bbc_ns_t[i] = std::numeric_limits<float>::quiet_NaN();
  }
}

void BbcVertexv2::identify(std::ostream& os) const
{
  os << "---BbcVertexv2--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << std::endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int BbcVertexv2::isValid() const
{
  if (_id == 0xFFFFFFFF)
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

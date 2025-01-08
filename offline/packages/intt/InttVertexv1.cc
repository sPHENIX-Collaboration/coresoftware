#include "InttVertexv1.h"

#include <cmath>

InttVertexv1::InttVertexv1()
  : _id(0xFFFFFFFF)
{
  std::fill(std::begin(_pos), std::end(_pos), std::numeric_limits<float>::quiet_NaN());
  std::fill(std::begin(_err), std::end(_err), std::numeric_limits<float>::quiet_NaN());
}

void InttVertexv1::identify(std::ostream& os) const
{
  os << "---InttVertexv1--------------------------------" << std::endl;
  os << "vertexid: " << get_id() << std::endl;

  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << std::endl;

  os << "         ( ";
  os << get_error(0, 0) << " , ";
  os << get_error(0, 1) << " , ";
  os << get_error(0, 2) << " )" << std::endl;
  os << "  err  = ( ";
  os << get_error(1, 0) << " , ";
  os << get_error(1, 1) << " , ";
  os << get_error(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << get_error(2, 0) << " , ";
  os << get_error(2, 1) << " , ";
  os << get_error(2, 2) << " )" << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int InttVertexv1::isValid() const
{
  if (_id == std::numeric_limits<unsigned int>::max())
  {
    return 0;
  }
  for (float _po : _pos)
  {
    if (!std::isfinite(_po))
    {
      return 0;
    }
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (!std::isfinite(get_error(i, j)))
      {
        return 0;
      }
    }
  }
  return 1;
}

void InttVertexv1::set_error(unsigned int i, unsigned int j, float value)
{
  _err[covar_index(i, j)] = value;
  return;
}

float InttVertexv1::get_error(unsigned int i, unsigned int j) const
{
  return _err[covar_index(i, j)];
}

unsigned int InttVertexv1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j)
  {
    std::swap(i, j);
  }
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

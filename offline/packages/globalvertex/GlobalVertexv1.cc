#include "GlobalVertexv1.h"

#include <cmath>

GlobalVertexv1::GlobalVertexv1(const GlobalVertex::VTXTYPE id)
  : _id(id)
{
  std::fill(std::begin(_pos), std::end(_pos), std::numeric_limits<float>::quiet_NaN());
  std::fill(std::begin(_err), std::end(_err), std::numeric_limits<float>::quiet_NaN());
}

void GlobalVertexv1::identify(std::ostream& os) const
{
  os << "---GlobalVertexv1-----------------------" << std::endl;
  os << "vertexid: " << get_id();

  os << " t = " << get_t() << std::endl;

  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << std::endl;

  os << " chisq = " << get_chisq() << ", ";
  os << " ndof = " << get_ndof() << std::endl;

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

  os << " list of vtx ids: " << std::endl;
  for (ConstVtxIter iter = begin_vtxids(); iter != end_vtxids(); ++iter)
  {
    os << "  " << iter->first << " => " << iter->second << std::endl;
  }
  os << "-----------------------------------------------" << std::endl;

  return;
}

int GlobalVertexv1::isValid() const
{
  if (!std::isfinite(_t))
  {
    return 0;
  }
  if (!std::isfinite(_t_err))
  {
    return 0;
  }
  if (!std::isfinite(_chisq))
  {
    return 0;
  }
  if (_ndof == std::numeric_limits<unsigned int>::max())
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
  if (_vtx_ids.empty())
  {
    return 0;
  }
  return 1;
}

void GlobalVertexv1::set_error(unsigned int i, unsigned int j, float value)
{
  _err[covar_index(i, j)] = value;
  return;
}

float GlobalVertexv1::get_error(unsigned int i, unsigned int j) const
{
  return _err[covar_index(i, j)];
}

unsigned int GlobalVertexv1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j)
  {
    std::swap(i, j);
  }
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

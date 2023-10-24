#include "GlobalVertexv2.h"

#include <cmath>

GlobalVertexv2::GlobalVertexv2(const GlobalVertex::VTXTYPE id)
  : _id(id)
{
  std::fill(std::begin(_pos), std::end(_pos), std::numeric_limits<float>::quiet_NaN());
  std::fill(std::begin(_err), std::end(_err), std::numeric_limits<float>::quiet_NaN());
}

void GlobalVertexv2::identify(std::ostream& os) const
{
  os << "---GlobalVertexv2-----------------------" << std::endl;
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

int GlobalVertexv2::isValid() const
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
  if (_vtxs.empty())
  {
    return 0;
  }
  return 1;
}
void GlobalVertexv2::insert_vtx(GlobalVertex::VTXTYPE type, Vertex* vertex)
{
  auto it = _vtxs.find(type);
  if( it == _vtxs.end())
    {
      VertexVector vector;
      vector.push_back(vertex);
      _vtxs.insert(std::make_pair(type, vector));
      return;
    }

  it->second.push_back(vertex);


}
size_t GlobalVertexv2::count_vtxs(GlobalVertex::VTXTYPE type) const
{
  auto it = _vtxs.find(type);
  if(it == _vtxs.end()) 
    {
      return 0;
    }

  return it->second.size();
}
void GlobalVertexv2::set_error(unsigned int i, unsigned int j, float value)
{
  _err[covar_index(i, j)] = value;
  return;
}

float GlobalVertexv2::get_error(unsigned int i, unsigned int j) const
{
  return _err[covar_index(i, j)];
}

unsigned int GlobalVertexv2::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j)
  {
    std::swap(i, j);
  }
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

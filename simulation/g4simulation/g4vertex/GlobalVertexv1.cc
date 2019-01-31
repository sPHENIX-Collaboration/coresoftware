#include "GlobalVertexv1.h"

#include <cmath>

using namespace std;

GlobalVertexv1::GlobalVertexv1()
  : _id(0xFFFFFFFF)
  , _t(NAN)
  , _t_err(NAN)
  , _pos()
  , _chisq(NAN)
  , _ndof(0xFFFFFFFF)
  , _err()
  , _vtx_ids()
{
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      set_error(i, j, NAN);
    }
  }
}

GlobalVertexv1::~GlobalVertexv1() {}

void GlobalVertexv1::identify(ostream& os) const
{
  os << "---GlobalVertexv1-----------------------" << endl;
  os << "vertexid: " << get_id() << endl;

  os << " t = " << get_t() << endl;

  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << endl;

  os << " chisq = " << get_chisq() << ", ";
  os << " ndof = " << get_ndof() << endl;

  os << "         ( ";
  os << get_error(0, 0) << " , ";
  os << get_error(0, 1) << " , ";
  os << get_error(0, 2) << " )" << endl;
  os << "  err  = ( ";
  os << get_error(1, 0) << " , ";
  os << get_error(1, 1) << " , ";
  os << get_error(1, 2) << " )" << endl;
  os << "         ( ";
  os << get_error(2, 0) << " , ";
  os << get_error(2, 1) << " , ";
  os << get_error(2, 2) << " )" << endl;

  os << " list of vtx ids: " << endl;
  for (ConstVtxIter iter = begin_vtxids(); iter != end_vtxids(); ++iter)
  {
    os << "  " << iter->first << " => " << iter->second << endl;
  }
  os << "-----------------------------------------------" << endl;

  return;
}

int GlobalVertexv1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  if (isnan(_t)) return 0;
  if (isnan(_t_err)) return 0;
  if (isnan(_chisq)) return 0;
  if (_ndof == 0xFFFFFFFF) return 0;

  for (int i = 0; i < 3; ++i)
  {
    if (isnan(_pos[i])) return 0;
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(get_error(i, j))) return 0;
    }
  }
  if (_vtx_ids.empty()) return 0;
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
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

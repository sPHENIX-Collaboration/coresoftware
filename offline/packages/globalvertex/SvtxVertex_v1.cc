#include "SvtxVertex_v1.h"

#include <cmath>
#include <utility>  // for swap

using namespace std;

SvtxVertex_v1::SvtxVertex_v1()
  : _id(0xFFFFFFFF)
  , _t0(NAN)
  , _pos()
  , _chisq(NAN)
  , _ndof(0xFFFFFFFF)
  , _err()
  , _track_ids()
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

void SvtxVertex_v1::identify(ostream& os) const
{
  os << "---SvtxVertex_v1--------------------" << endl;
  os << "vertexid: " << get_id() << endl;

  os << " t0 = " << get_t0() << endl;

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

  os << " list of tracks ids: ";
  for (ConstTrackIter iter = begin_tracks(); iter != end_tracks(); ++iter)
  {
    os << *iter << " ";
  }
  os << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int SvtxVertex_v1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  if (std::isnan(_t0)) return 0;
  if (std::isnan(_chisq)) return 0;
  if (_ndof == 0xFFFFFFFF) return 0;

  for (int i = 0; i < 3; ++i)
  {
    if (std::isnan(_pos[i])) return 0;
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (std::isnan(get_error(i, j))) return 0;
    }
  }
  if (_track_ids.empty()) return 0;
  return 1;
}

void SvtxVertex_v1::set_error(unsigned int i, unsigned int j, float value)
{
  _err[covar_index(i, j)] = value;
  return;
}

float SvtxVertex_v1::get_error(unsigned int i, unsigned int j) const
{
  return _err[covar_index(i, j)];
}

unsigned int SvtxVertex_v1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

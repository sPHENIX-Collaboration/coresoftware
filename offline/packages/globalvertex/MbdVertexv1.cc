#include "MbdVertexv1.h"

#include <cmath>

using namespace std;

MbdVertexv1::MbdVertexv1()
  : _id(0xFFFFFFFF)
  , _t(NAN)
  , _t_err(NAN)
  , _z(NAN)
  , _z_err(NAN)
{
}

MbdVertexv1::~MbdVertexv1() = default;

void MbdVertexv1::identify(ostream& os) const
{
  os << "---MbdVertexv1--------------------------------" << endl;
  os << "vertexid: " << get_id() << endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int MbdVertexv1::isValid() const
{
  if (_id == 0xFFFFFFFF)
  {
    return 0;
  }
  if (isnan(_t))
  {
    return 0;
  }
  if (isnan(_t_err))
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

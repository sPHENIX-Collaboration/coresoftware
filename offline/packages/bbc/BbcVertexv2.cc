#include "BbcVertexv2.h"

#include <cmath>

using namespace std;

BbcVertexv2::BbcVertexv2()
  : _id(0xFFFFFFFF)
  , _t(NAN)
  , _t_err(NAN)
  , _z(NAN)
  , _z_err(NAN)
{
  for (int i = 0; i < 2; i++)
  {
    _bbc_ns_npmt[i] = 0;
    _bbc_ns_q[i] = NAN;
    _bbc_ns_t[i] = NAN;
  }
}

BbcVertexv2::~BbcVertexv2() = default;

void BbcVertexv2::identify(ostream& os) const
{
  os << "---BbcVertexv2--------------------------------" << endl;
  os << "vertexid: " << get_id() << endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int BbcVertexv2::isValid() const
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

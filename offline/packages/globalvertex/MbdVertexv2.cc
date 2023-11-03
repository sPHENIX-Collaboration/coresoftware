#include "MbdVertexv2.h"

#include <cmath>

using namespace std;

MbdVertexv2::MbdVertexv2()
  : _id(0xFFFFFFFF)
  , _bco(std::numeric_limits<unsigned int>::max())
  , _t(NAN)
  , _t_err(NAN)
  , _z(NAN)
  , _z_err(NAN)
{
}

MbdVertexv2::~MbdVertexv2() = default;

void MbdVertexv2::identify(ostream& os) const
{
  os << "---MbdVertexv2--------------------------------" << endl;
  os << "vertexid: " << get_id() << endl;
  os << " t = " << get_t() << " +/- " << get_t_err() << endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int MbdVertexv2::isValid() const
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

float MbdVertexv2::get_position(unsigned int coor) const 
{
  if(coor == 0) return get_x();
  else if(coor == 1) return get_y();
  else if(coor == 2) return get_z();
  else return NAN;
}

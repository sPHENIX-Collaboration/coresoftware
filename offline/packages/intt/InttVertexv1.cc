#include "InttVertexv1.h"

#include <cmath>

using namespace std;

InttVertexv1::InttVertexv1()
  : _id(0xFFFFFFFF)
  , _z(NAN)
  , _z_err(NAN)
{
}

InttVertexv1::~InttVertexv1() = default;

void InttVertexv1::identify(ostream& os) const
{
  os << "---InttVertexv1--------------------------------" << endl;
  os << "vertexid: " << get_id() << endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int InttVertexv1::isValid() const
{
  if (_id == 0xFFFFFFFF)
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

#include "SimpleHit3D.h"

#include <iostream>

using namespace std;

SimpleHit3D::SimpleHit3D(float xx, float dxx,
			 float yy, float dyy,
			 float zz, float dzz,
			 unsigned int ind, int lyr) 
  : x(xx), dx(dxx),
    y(yy), dy(dyy),
    z(zz), dz(dzz),
    index(ind), layer(lyr),
    _err()
{
  for (int i = 0; i < 3; ++i) _err[i] = new float[i+1];
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,0.0);
    }
  }
}			 

SimpleHit3D::SimpleHit3D(const SimpleHit3D &hit)
  : x(hit.x), dx(hit.dx),
    y(hit.y), dy(hit.dy),
    z(hit.z), dz(hit.dz),
    index(hit.index), layer(hit.layer),
    _err()
{
  for (int i=0; i<3; ++i) _err[i] = new float[i+1];
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_error(i,j,hit.get_error(i,j));
    }
  } 
}

SimpleHit3D& SimpleHit3D::operator= (const SimpleHit3D &rhs) {
  x = rhs.x; dx = rhs.dx;
  y = rhs.y; dy = rhs.dy;
  z = rhs.z; dz = rhs.dz;
  index = rhs.index;
  layer = rhs.layer;
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_error(i,j,rhs.get_error(i,j));
    }
  } 
  
  return *this;
}

SimpleHit3D::~SimpleHit3D() {
  for (int i=0; i<3; ++i) delete[] _err[i];
}

void SimpleHit3D::print(std::ostream& out) const {

  out << "SimpleHit3D: "
       << "id: " << index << " layer: " << layer << " "
       << "(x,y,z) = (" << x << "," << y << "," << z << ") "
       << "(dx,dy,dz) = (" << dx << "," << dy << "," << dz << ")" << endl;

  out << "       ( ";
  out << get_error(0,0) << " , ";
  out << get_error(0,1) << " , ";
  out << get_error(0,2) << " )" << endl; 
  out << "err  = ( ";
  out << get_error(1,0) << " , ";
  out << get_error(1,1) << " , ";
  out << get_error(1,2) << " )" << endl;
  out << "       ( ";
  out << get_error(2,0) << " , ";
  out << get_error(2,1) << " , ";
  out << get_error(2,2) << " )" << endl;

  return;
}

void SimpleHit3D::set_error(int i, int j, float value) {
  if (j > i) set_error(j,i,value);
  else _err[i][j] = value;
  return;
}

float SimpleHit3D::get_error(int i, int j) const {
  if (j > i) return get_error(j,i);
  return _err[i][j];
}

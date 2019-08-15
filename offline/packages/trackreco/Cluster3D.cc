#include "Cluster3D.h"

#include <iostream>
#include <utility>   // for swap

using namespace std;


Cluster3D::Cluster3D()
  : _cluskey(0),
    _id(0),
    _layer(-1),
    _x(0.0),
    _y(0.0),
    _z(0.0),
    _err(),
    _size()
{
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,0.0);
      set_size(i,j,0.0);
    }
  }
}	

void Cluster3D::print(std::ostream& out) const {

  out << "Cluster3D: "
      << "id: " << get_id() << " cluskey " << _cluskey << " layer: " << get_layer() << " "
      << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ") "
      << endl;

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

  out << "       ( ";
  out << get_size(0,0) << " , ";
  out << get_size(0,1) << " , ";
  out << get_size(0,2) << " )" << endl; 
  out << "size = ( ";
  out << get_size(1,0) << " , ";
  out << get_size(1,1) << " , ";
  out << get_size(1,2) << " )" << endl;
  out << "       ( ";
  out << get_size(2,0) << " , ";
  out << get_size(2,1) << " , ";
  out << get_size(2,2) << " )" << endl;
  
  return;
}

void Cluster3D::set_error(unsigned int i, unsigned int j, float value) {
  _err[covar_index(i,j)] = value;
  return;
}

float Cluster3D::get_error(unsigned int i, unsigned int j) const {
  return _err[covar_index(i,j)];
}

void Cluster3D::set_size(unsigned int i, unsigned int j, float value) {
  _size[covar_index(i,j)] = value;
  return;
}

float Cluster3D::get_size(unsigned int i, unsigned int j) const {
  return _size[covar_index(i,j)];
}

unsigned int Cluster3D::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}

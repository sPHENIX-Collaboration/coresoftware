#ifndef __CLUSTER3D_H__
#define __CLUSTER3D_H__

#include <iostream>
#include <cmath>

class Cluster3D
{

public:

  Cluster3D();
  virtual ~Cluster3D() {}

  unsigned int get_id() const {return _id;}
  void         set_id(unsigned int id) {_id = id;}

  int  get_layer() const {return _layer;}
  void set_layer(int layer) {_layer = layer;}
  
  float get_x() const {return _x;}
  void  set_x(float x) {_x = x;}

  float get_y() const {return _y;}
  void  set_y(float y) {_y = y;}

  float get_z() const {return _z;}
  void  set_z(float z) {_z = z;}
  
  void  print(std::ostream& out = std::cout) const; //< dump the values to screen

  float get_error(unsigned int i, unsigned int j) const;
  void  set_error(unsigned int i, unsigned int j, float value);

  float get_size(unsigned int i, unsigned int j) const;
  void  set_size(unsigned int i, unsigned int j, float value);
   
private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id;
  int _layer;
  
  float _x;
  float _y;
  float _z;
  
  float _err[6]; //< error covariance matrix (x,y,z)
  float _size[6]; //< size covariance matrix (x,y,z)
};

#endif // __CLUSTER3D_H__

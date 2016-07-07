#ifndef __SIMPLEHIT3D__
#define __SIMPLEHIT3D__

#include <iostream>

class SimpleHit3D
{

public:
  
  SimpleHit3D(float x = 0.0, float dxx = 0.0,
	      float y = 0.0, float dyy = 0.0,
	      float z = 0.0, float dzz = 0.0,
	      unsigned int id = 0, int lyr = -1);
  virtual ~SimpleHit3D() {}

  unsigned int get_id() const {return _id;}
  void         set_id(unsigned int id) {_id = id;}

  float get_x() const {return _x;}
  void  set_x(float x) {_x = x;}

  float get_y() const {return _y;}
  void  set_y(float y) {_y = y;}

  float get_z() const {return _z;}
  void  set_z(float z) {_z = z;}
  
  void  print(std::ostream& out = std::cout) const; //< dump the values to screen

  float get_error(unsigned int i, unsigned int j) const;        //< get cluster error covar
  void  set_error(unsigned int i, unsigned int j, float value); //< set cluster error covar
  
private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id;

  float _x;
  float _y;
  float _z;
  
public:
  float dx;
  float dy;
  float dz;  
  int layer;

private:
  
  float _err[6]; //< error covariance matrix (x,y,z) (jagged array) 
};

#endif // __SIMPLEHIT3D__

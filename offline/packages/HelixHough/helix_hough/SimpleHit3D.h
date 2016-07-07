#ifndef __SIMPLEHIT3D__
#define __SIMPLEHIT3D__

#include <iostream>

class SimpleHit3D
{

public:
  
  SimpleHit3D(float xx = 0.0, float dxx = 0.0,
	      float yy = 0.0, float dyy = 0.0,
	      float zz = 0.0, float dzz = 0.0,
	      unsigned int id = 0, int lyr = -1);
  virtual ~SimpleHit3D() {}

  unsigned int get_id() const {return _id;}
  void         set_id(unsigned int id) {_id = id;}

  void  print(std::ostream& out = std::cout) const; //< dump the values to screen

  float get_size(unsigned int i, unsigned int j) const;        //< get cluster size covar
  void  set_size(unsigned int i, unsigned int j, float value); //< set cluster size covar
  
  float get_error(unsigned int i, unsigned int j) const;        //< get cluster error covar
  void  set_error(unsigned int i, unsigned int j, float value); //< set cluster error covar
  
private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id;

public:
  float x, dx;
  float y, dy;
  float z, dz;  
  int layer;

private:
  
  float _size[6]; //< size covariance matrix (x,y,z) (jagged array) 
  float _err[6]; //< error covariance matrix (x,y,z) (jagged array) 
};

#endif // __SIMPLEHIT3D__

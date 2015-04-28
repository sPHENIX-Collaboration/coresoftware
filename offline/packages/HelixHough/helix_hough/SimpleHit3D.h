#ifndef __SIMPLEHIT3D__
#define __SIMPLEHIT3D__

#include <iostream>

class SimpleHit3D
{

public:
  
  SimpleHit3D(float xx=0.,float dxx=0.,
	      float yy=0., float dyy=0.,
	      float zz=0., float dzz=0.,
	      unsigned int ind=0, int lyr=-1);
  SimpleHit3D(const SimpleHit3D &hit);
  SimpleHit3D& operator=(const SimpleHit3D &hit);
  virtual ~SimpleHit3D();

  float x, dx;
  float y, dy;
  float z, dz;  
  unsigned int index;
  int layer;

  void  print(std::ostream& out = std::cout) const; //< dump the values to screen
  
  float get_error(int i, int j) const;        //< get cluster error covar
  void  set_error(int i, int j, float value); //< set cluster error covar
  
private:

  float* _err[3]; //< error covariance matrix (x,y,z) (jagged array) 
};

#endif // __SIMPLEHIT3D__

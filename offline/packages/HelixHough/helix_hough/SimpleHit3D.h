#ifndef __SIMPLEHIT3D__
#define __SIMPLEHIT3D__

#include <iostream>
#include <cmath>
#include <trackbase/TrkrDefs.h>

class SimpleHit3D
{

public:

  SimpleHit3D();
  virtual ~SimpleHit3D() {}

  TrkrDefs::cluskey get_cluskey() const {return _cluskey;}
  void         set_cluskey(TrkrDefs::cluskey key) {_cluskey = key;}

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

  TrkrDefs::cluskey _cluskey;
  unsigned int _id;
  int _layer;
  
  float _x;
  float _y;
  float _z;
  
  float _err[6]; //< error covariance matrix (x,y,z)
  float _size[6]; //< size covariance matrix (x,y,z)
  
};

#endif // __SIMPLEHIT3D__

#ifndef __SVTXBEAMSPOT_H__
#define __SVTXBEAMSPOT_H__

#include <phool/PHObject.h>
#include <vector>
#include <set>
#include <iostream>

class SvtxBeamSpot : public PHObject {

public:
  
  SvtxBeamSpot();
  SvtxBeamSpot(const SvtxBeamSpot &beamspot);
  SvtxBeamSpot& operator=(const SvtxBeamSpot &beamspot);
  virtual ~SvtxBeamSpot();

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          IsValid() const;

  // beamspot info
 
  float        get_x() const                         {return _pos[0];}
  void         set_x(float x)                        {_pos[0] = x;}
  
  float        get_y() const                         {return _pos[1];}
  void         set_y(float y)                        {_pos[1] = y;}
  
  float        get_position(int coor) const          {return _pos[coor];}
  void         set_position(int coor, float xi)      {_pos[coor] = xi;}

  float        get_error(int i, int j) const;        //< get beamspot error covar
  void         set_error(int i, int j, float value); //< set beamspot error covar
  
private:

  float  _pos[2]; //< position x,y
  float* _err[2]; //< variance covariance matrix (+/- cm^2)
  
  ClassDef(SvtxBeamSpot, 1);
};

#endif


#ifndef G4HOUGH_SVTXBEAMSPOT_H
#define G4HOUGH_SVTXBEAMSPOT_H

#include <phool/PHObject.h>

#include <iostream>

class SvtxBeamSpot : public PHObject {

public:
  
  SvtxBeamSpot();
  virtual ~SvtxBeamSpot() {}

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset() {*this = SvtxBeamSpot();}
  int          isValid() const;

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

  unsigned int covar_index(unsigned int i, unsigned int j) const;
  
  float  _pos[2]; //< position x,y
  float  _err[3]; //< variance covariance matrix (packed storage) (+/- cm^2)
  
  ClassDef(SvtxBeamSpot, 1);
};

#endif


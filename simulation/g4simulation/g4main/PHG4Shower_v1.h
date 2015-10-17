#ifndef __PHG4SHOWER_V1_H__
#define __PHG4SHOWER_V1_H__

#include "PHG4Shower.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class PHG4Shower_v1 : public PHG4Shower {

public:
  
  PHG4Shower_v1();
  virtual ~PHG4Shower_v1();

  // PHObject virtual overloads
   
  void         identify(std::ostream& os = std::cout) const;
  PHG4Shower*  Clone();
  void         Reset();
  int          isValid() const;

  // shower info
  
  unsigned int get_id() const                 {return _id;}
  void         set_id(unsigned int id)        {_id = id;}

  int          get_primary_id() const         {return _primary_id;}
  void         set_primary_id(int primary_id) {_primary_id = primary_id;}
  
  float        get_x() const                  {return _pos[0];}
  void         set_x(float x)                 {_pos[0] = x;}

  float        get_y() const                  {return _pos[1];}
  void         set_y(float y)                 {_pos[1] = y;}

  float        get_z() const                  {return _pos[2];}
  void         set_z(float z)                 {_pos[2] = z;}

  float        get_position(unsigned int coor) const     {return _pos[coor];}
  void         set_position(unsigned int coor, float xi) {_pos[coor] = xi;}
  
  float        get_covar(unsigned int i, unsigned int j) const;
  void         set_covar(unsigned int i, unsigned int j, float entry);
  
  float        get_edep(PHG4Shower::VOLUME calotype) const;
  void         set_edep(PHG4Shower::VOLUME calotype, float edep);

  float        get_eion(PHG4Shower::VOLUME calotype) const;
  void         set_eion(PHG4Shower::VOLUME calotype, float eion);

  float        get_light_yield(PHG4Shower::VOLUME calotype);
  void         set_light_yield(PHG4Shower::VOLUME calotype, float light_yield);
  
private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;
  
  unsigned int                        _id;          //< unique identifier within container
  float                               _pos[3];      //< mean position of the shower
  float                               _covar[6];    //< covariance of shower positions
  std::map<PHG4Shower::VOLUME, float> _edep;        //< energy deposit in different volumes
  std::map<PHG4Shower::VOLUME, float> _eion;        //< ionization energy in different volumes
  std::map<PHG4Shower::VOLUME, float> _light_yield; //< light yield in different volumes

  
  ClassDef(PHG4Shower_v1, 1);
};

#endif


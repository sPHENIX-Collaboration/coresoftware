#ifndef __PHG4SHOWER_H__
#define __PHG4SHOWER_H__

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>
#include <cmath>
#include <set>
#include <iostream>
#include <string>

class PHG4Shower : public PHObject {

public:

  virtual ~PHG4Shower() {}

  // PHObject virtual overloads
  
  virtual void         identify(std::ostream& os = std::cout) const {os << "PHG4Shower base class" << std::endl;}
  virtual PHG4Shower*  Clone() const                                {return NULL;}
  virtual void         Reset()                                      {}
  virtual int          isValid() const                              {return 0;}

  // shower info
  
  virtual unsigned int get_id() const           {return 0xFFFFFFFF;}
  virtual void         set_id(unsigned int id)  {}

  virtual int          get_primary_id() const           {return -1;}
  virtual void         set_primary_id(int primary_id)   {}
  
  virtual float        get_x() const            {return NAN;}
  virtual void         set_x(float x)           {}

  virtual float        get_y() const            {return NAN;}
  virtual void         set_y(float y)           {}

  virtual float        get_z() const            {return NAN;}
  virtual void         set_z(float x)           {}

  virtual float        get_position(unsigned int coor) const          {return NAN;}
  virtual void         set_position(unsigned int coor, float xi)      {}
  
  virtual float        get_covar(unsigned int i, unsigned int j) const {return NAN;}
  virtual void         set_covar(unsigned int i, unsigned int j, float entry) {}
  
  virtual float        get_edep(int volume) const {return NAN;}
  virtual void         set_edep(int volume, float edep) {}

  virtual float        get_eion(int volume) const {return NAN;}
  virtual void         set_eion(int volume, float eion) {}

  virtual float        get_light_yield(int volume) const {return NAN;}
  virtual void         set_light_yield(int volume, float light_yield) {}

  virtual void         add_g4particle_id(int id)  {}
  virtual size_t       remove_g4particle_id(int id) {return 0;}
  
  virtual void         add_g4hit_id(int volume, PHG4HitDefs::keytype id) {}
  virtual size_t       remove_g4hit_id(int volume, PHG4HitDefs::keytype id) {return 0;}
  
protected:
  PHG4Shower() {}
  
private:
  
  ClassDef(PHG4Shower, 1);
};

#endif


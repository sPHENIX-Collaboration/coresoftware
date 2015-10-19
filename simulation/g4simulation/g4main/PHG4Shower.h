#ifndef __PHG4SHOWER_H__
#define __PHG4SHOWER_H__

#include <phool/PHObject.h>
#include <cmath>
#include <iostream>

class PHG4Shower : public PHObject {

public:

  enum VOLUME {NONE=0,
	       CEMC_ELECTRONICS=1,CEMC=2,ABSORBER_CEMC=3,CEMC_SPT=4,
	       ABSORBER_HCALIN=5,HCALIN=6,HCALIN_SPT=7,
	       MAGNET=8,
	       ABSORBER_HCALOUT=9,HCALOUT=10,HCALOUT_SPT=11,
	       BH_1=12
  };
  
  virtual ~PHG4Shower() {}

  // PHObject virtual overloads
  
  virtual void         identify(std::ostream& os = std::cout) const {os << "PHG4Shower base class" << std::endl;}
  virtual PHG4Shower*  Clone()                                      {return NULL;}
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
  
  virtual float        get_edep(VOLUME calotype) const {return NAN;}
  virtual void         set_edep(VOLUME calotype, float edep) {}

  virtual float        get_eion(VOLUME calotype) const {return NAN;}
  virtual void         set_eion(VOLUME calotype, float eion) {}

  virtual float        get_light_yield(VOLUME calotype) const {return NAN;}
  virtual void         set_light_yield(VOLUME calotype, float light_yield) {}
  
protected:
  PHG4Shower() {}
  
private:
  
  ClassDef(PHG4Shower, 1);
};

#endif


/*
 Originated by Tristan Protzman 12/15/22
 Re-written by Ejiro Umaka 03/28/23
*/

#ifndef EPD_EPDGEOM_H
#define EPD_EPDGEOM_H

#include <phool/PHObject.h>

#include <limits>

class EpdGeom : public PHObject
{
  public:
    EpdGeom() = default;
    ~EpdGeom() override {};

    virtual void set_z(unsigned int /*key*/, float /*z*/) {return;}
    virtual void set_r(unsigned int /*key*/, float /*r*/) {return;}
    virtual void set_phi(unsigned int /*key*/, float /*f*/) {return;}
    virtual void set_phi0(unsigned int /*key*/, float /*f0*/) {return;} 
    virtual float get_r(unsigned int /*key*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_z(unsigned int /*key*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_phi(unsigned int /*key*/) const {return std::numeric_limits<float>::quiet_NaN();};

  private:
    ClassDefOverride(EpdGeom, 1);
};


#endif // EPD_EPDGEOM_H

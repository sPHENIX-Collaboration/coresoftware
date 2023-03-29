/*
 Originated by Tristan Protzman 12/15/22
 Re-written by Ejiro Umaka 03/28/23
*/

#ifndef EPD_GEOM_H
#define EPD_GEOM_H

#include <vector>
#include <utility>
#include <tuple>
#include <iostream>

#include <phool/PHObject.h>

class EpdGeom : public PHObject
{
  public:
    EpdGeom() = default;
    ~EpdGeom() override {};

    virtual unsigned int get_arm_index(unsigned int /*key*/){return 999;};
    virtual unsigned int get_r_index(unsigned int /*key*/){return 999;};
    virtual unsigned int get_phi_index(unsigned int /*key*/){return 999;};
    virtual void set_z(unsigned int /*key*/, float /*z*/) {return;}
    virtual void set_r(unsigned int /*key*/, float /*r*/) {return;}
    virtual void set_phi(unsigned int /*key*/, float /*f*/) {return;}
    virtual float get_r(unsigned int /*key*/){return 999;};
    virtual float get_z(unsigned int /*key*/){return 999;};
    virtual float get_phi(unsigned int /*key*/){return 999;};

    virtual void Reset() override {return;}; 

  private:
    ClassDefOverride(EpdGeom, 1);
};


#endif // EPD_GEOM_H

#ifndef __BBC_GEOM_H__
#define __BBC_GEOM_H__

#include <phool/PHObject.h>

#include <vector>
#include <utility>
#include <tuple>
#include <iostream>
#include <cmath>

class BbcGeom : public PHObject
{
  public:
    BbcGeom() = default;
    ~BbcGeom() override {};

    virtual float get_x(const unsigned int /*pmtch*/) const {return NAN;};
    virtual float get_y(const unsigned int /*pmtch*/) const {return NAN;};
    virtual float get_z(const unsigned int /*pmtch*/) const {return NAN;};
    virtual float get_r(const unsigned int /*pmtch*/) const {return NAN;};
    virtual float get_phi(const unsigned int /*pmtch*/) const {return NAN;};
    virtual int   get_arm(const unsigned int /*pmtch*/) const {return -1;};
    virtual void  set_xyz(const unsigned int /*pmtch*/, const float /*x*/, const float /*y*/, const float /*z*/) {}

    // methods when accessing raw fee channels
    virtual int get_arm_feech(const unsigned int /*feech*/) const {return -1;};
    virtual int get_pmt(const unsigned int /*feech*/) const {return -1;};
    virtual int get_type(const unsigned int /*feech*/) const {return -1;}; // 0=T-channel, 1=Q-channel

  private:
    ClassDefOverride(BbcGeom, 1);
};


#endif // __BBC_GEOM_H__

#ifndef MBD_GEOM_H
#define MBD_GEOM_H

#include <phool/PHObject.h>

#include <vector>
#include <utility>
#include <tuple>
#include <iostream>
#include <limits>

class MbdGeom : public PHObject
{
  public:
    MbdGeom() = default;
    ~MbdGeom() override = default;

    virtual float get_x(const unsigned int /*pmtch*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_y(const unsigned int /*pmtch*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_z(const unsigned int /*pmtch*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_r(const unsigned int /*pmtch*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual float get_phi(const unsigned int /*pmtch*/) const {return std::numeric_limits<float>::quiet_NaN();};
    virtual int   get_arm(const unsigned int pmtch) const { return pmtch / 64; }
    virtual int   get_feech(const unsigned int pmtch, const unsigned int type = 1) const {
      return (pmtch / 8) * 16 + pmtch % 8 + 8*type;
    }
    virtual void  set_xyz(const unsigned int /*pmtch*/, const float /*x*/, const float /*y*/, const float /*z*/) {}

    // methods when accessing raw fee channels
    virtual int get_arm_feech(const unsigned int feech) const { return get_pmt(feech) / 64; }
    virtual int get_pmt(const unsigned int feech) const { return (feech / 16) * 8 + feech % 8; }
    virtual int get_type(const unsigned int feech) const { return (feech / 8) % 2; }  // 0=T-channel, 1=Q-channel

    virtual void download_hv() {}

    virtual void Reset() override {}

  private:
    ClassDefOverride(MbdGeom, 1);
};


#endif // __MBD_GEOM_H__

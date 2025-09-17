#ifndef CALOBASE_RAWCLUSTERV2_H
#define CALOBASE_RAWCLUSTERV2_H

#include "RawClusterv1.h"
#include <limits>

class RawClusterv2 : public RawClusterv1
{
 public:
  RawClusterv2() = default;
  ~RawClusterv2() override = default;

  // ---- new API: tower-space CoG, in TOWER UNITS ----
  void set_tower_cog(float xr, float yr, float xc, float yc)
  { _x_raw = xr; _y_raw = yr; _x_corr = xc; _y_corr = yc; }

  float x_tower_raw()  const { return _x_raw;  }
  float y_tower_raw()  const { return _y_raw;  }
  float x_tower_corr() const { return _x_corr; }
  float y_tower_corr() const { return _y_corr; }

  // Optional: clone into the same dynamic type
  RawCluster* CloneMe() const override { return new RawClusterv2(*this); }

  void Reset() override
  {
    RawClusterv1::Reset();
    _x_raw  = _y_raw  = _x_corr = _y_corr = std::numeric_limits<float>::quiet_NaN();
  }

  void identify(std::ostream& os = std::cout) const override
  {
    RawClusterv1::identify(os);
    os << "  [towerCoG raw=(" << _x_raw << "," << _y_raw << ") corr=("
       << _x_corr << "," << _y_corr << ")]\n";
  }

 private:
  float _x_raw  = std::numeric_limits<float>::quiet_NaN();
  float _y_raw  = std::numeric_limits<float>::quiet_NaN();
  float _x_corr = std::numeric_limits<float>::quiet_NaN();
  float _y_corr = std::numeric_limits<float>::quiet_NaN();

  ClassDefOverride(RawClusterv2, 1)
};

#endif


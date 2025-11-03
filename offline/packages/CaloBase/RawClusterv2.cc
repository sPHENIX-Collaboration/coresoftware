#include "RawClusterv2.h"

void RawClusterv2::Reset()
{
  RawClusterv1::Reset();
  _x_raw = std::numeric_limits<float>::quiet_NaN();
  _y_raw = std::numeric_limits<float>::quiet_NaN();
  _x_corr = std::numeric_limits<float>::quiet_NaN();
  _y_corr = std::numeric_limits<float>::quiet_NaN();
  _t_mean = std::numeric_limits<float>::quiet_NaN();
}

void RawClusterv2::identify(std::ostream& os) const
{
  RawClusterv1::identify(os);
  os << "  [towerCoG raw=(" << _x_raw << "," << _y_raw << ") corr=("
     << _x_corr << "," << _y_corr << ")]";
  os << "  [t_mean=" << _t_mean << "]\n";
}

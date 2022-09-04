#include "MapToPadPlanePassData.h"

class MapToPadPlanePassData;

MapToPadPlanePassData& MapToPadPlanePassData::operator+=(const MapToPadPlanePassData& rhs) {
  hitsetkey = rhs.hitsetkey;
  if (rhs.phi_bin_lo < phi_bin_lo) phi_bin_lo = rhs.phi_bin_lo;
  if (rhs.phi_bin_hi > phi_bin_hi) phi_bin_hi = rhs.phi_bin_hi;

  if (rhs.time_bin_lo < time_bin_lo) time_bin_lo = rhs.time_bin_lo;
  if (rhs.time_bin_hi > time_bin_hi) time_bin_hi = rhs.time_bin_hi;

  neff_electrons += rhs.neff_electrons;
  phi_integral   += rhs.phi_integral;
  time_integral  += rhs.time_integral;

  return *this;
}

void MapToPadPlanePassData::reset() {
  hitsetkey = 0;
  neff_electrons = 0;
  phi_integral = 0.;
  time_integral = 0.;
  phi_bin_lo = 0;
  phi_bin_hi = 0;
  time_bin_lo = 0;
  time_bin_hi = 0;
}

MapToPadPlanePassData::operator short() {
  return TrkrDefs::getLayer(hitsetkey);
}

bool MapToPadPlanePassData::has_data() const {
  return neff_electrons != 0;
}


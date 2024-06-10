#ifndef __MBD_GEOM_V2_H__
#define __MBD_GEOM_V2_H__

#include "MbdGeom.h"
#include <map>

class MbdGeomV2 : public MbdGeom
{
 public:
  MbdGeomV2();
  ~MbdGeomV2() override = default;

  float get_x(const unsigned int pmtch) const override { return pmt_x[pmtch]; }
  float get_y(const unsigned int pmtch) const override { return pmt_y[pmtch]; }
  float get_z(const unsigned int pmtch) const override { return pmt_z[pmtch]; }
  float get_r(const unsigned int pmtch) const override { return pmt_r[pmtch]; }
  float get_phi(const unsigned int pmtch) const override { return pmt_phi[pmtch]; }

  void set_xyz(const unsigned int ipmt, const float x, const float y, const float z) override;

  void download_hv() override;
  const std::multimap<int,int>& get_hvmap();

  virtual void Reset() override {}

 private:
  float pmt_x[128]{};
  float pmt_y[128]{};
  float pmt_z[128]{};
  float pmt_r[128]{};
  float pmt_phi[128]{};

  std::multimap<int,int> pmt_hv;

  ClassDefOverride(MbdGeomV2, 1)
};

#endif  // __MBD_GEOM_V2_H__

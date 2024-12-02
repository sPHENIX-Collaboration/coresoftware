#ifndef __MBD_GEOM_V1_H__
#define __MBD_GEOM_V1_H__

#include "MbdGeom.h"

class MbdGeomV1 : public MbdGeom
{
 public:
  MbdGeomV1();
  ~MbdGeomV1() override = default;

  float get_x(const unsigned int pmtch) const override { return pmt_x[pmtch]; }
  float get_y(const unsigned int pmtch) const override { return pmt_y[pmtch]; }
  float get_z(const unsigned int pmtch) const override { return pmt_z[pmtch]; }
  float get_r(const unsigned int pmtch) const override { return pmt_r[pmtch]; }
  float get_phi(const unsigned int pmtch) const override { return pmt_phi[pmtch]; }

  void set_xyz(const unsigned int ipmt, const float x, const float y, const float z) override;

  virtual void Reset() override {}

 private:
  float pmt_x[128]{};
  float pmt_y[128]{};
  float pmt_z[128]{};
  float pmt_r[128]{};
  float pmt_phi[128]{};

  ClassDefOverride(MbdGeomV1, 1)
};

#endif  // __MBD_GEOM_V1_H__

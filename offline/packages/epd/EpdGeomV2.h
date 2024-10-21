#ifndef EPD_EPDGEOMV2_H
#define EPD_EPDGEOMV2_H

#include "EpdGeom.h"

#include <array>

class EpdGeomV2 : public EpdGeom {

public:
  EpdGeomV2() = default;
  ~EpdGeomV2() override = default;
  
  float get_r(unsigned int key) const override;
  float get_z(unsigned int key) const override;
  float get_phi(unsigned int key) const override;
  void set_z(unsigned int key, float z) override;
  void set_r(unsigned int key, float r) override;
  void set_phi(unsigned int key, float f) override;
  void set_phi0(unsigned int key, float f0) override;


private:
  
  std::array<float,16> tile_r {};
  std::array<float,2> tile_z {};
  std::array<float,24> tile_phi {};
  std::array<float,12> tile_phi0 {};

  ClassDefOverride(EpdGeomV2, 1)

};

#endif // EPD_EPDGEOMV2_H

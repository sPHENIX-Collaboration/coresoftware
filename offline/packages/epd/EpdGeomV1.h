#ifndef EPD_GEOM_V1_H
#define EPD_GEOM_V1_H

#include "EpdGeom.h"

#include <utility>
#include <tuple>
#include <iostream>


class EpdGeomV1 : public EpdGeom {

public:
  EpdGeomV1() = default;
  ~EpdGeomV1() override = default;
  
  float get_r(unsigned int key) const override;
  float get_z(unsigned int key) const override;
  float get_phi(unsigned int key) const override;
  void set_z(unsigned int key, float z) override;
  void set_r(unsigned int key, float r) override;
  void set_phi(unsigned int key, float f) override;
  void set_phi0(unsigned int key, float f0) override;


private:
  
  float tile_r[16] = {};
  float tile_z[2] = {};
  float tile_phi[24] = {};
  float tile_phi0[12] = {};

  ClassDefOverride(EpdGeomV1, 1)

};

#endif // EPD_GEOM_V1_H

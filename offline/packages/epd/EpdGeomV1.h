#ifndef EPD_GEOM_V1_H
#define EPD_GEOM_V1_H

#include "EpdGeom.h"
#include "EpdGeomV1.h"

#include <map>
#include <utility>
#include <tuple>
#include <iostream>

#include <phool/PHObject.h>


class EpdGeomV1 : public EpdGeom {

public:
  EpdGeomV1();
  ~EpdGeomV1() override;
  
  unsigned int get_arm_index(unsigned int key) override;
  unsigned int get_r_index(unsigned int key) override;
  unsigned int get_phi_index(unsigned int key) override;
  float get_r(unsigned int key) override {return tile_r[get_r_index(key)];}
  float get_z(unsigned int key) override {return tile_z[get_arm_index(key)];}
  float get_phi(unsigned int key) override {return tile_phi[key];}
  void set_z(unsigned int key, float z) override {tile_z[get_arm_index(key)] = z;}
  void set_r(unsigned int key, float r) override {tile_r[get_r_index(key)] = r;}
  void set_phi(unsigned int key, float f) override {tile_phi[key] = f;}

  void Reset() override {return;}; 

private:
  
  float tile_r[16] = {};
  float tile_z[2] = {};
  std::map <int, float> tile_phi;
  
  ClassDefOverride(EpdGeomV1, 1)

};

#endif // EPD_GEOM_V1_H

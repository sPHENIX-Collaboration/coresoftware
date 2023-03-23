#ifndef EPD_GEOM_H
#define EPD_GEOM_H

#include <vector>
#include <utility>
#include <tuple>
#include <iostream>

#include <phool/PHObject.h>

class EpdGeom : public PHObject
{
  public:
    EpdGeom() = default;
    ~EpdGeom() override {};

    virtual unsigned int side_r_phi_to_id(unsigned int /*side*/, unsigned int /*r_index*/, unsigned int /*phi_index*/) {return 999;};
    virtual unsigned int side_sector_tile_to_id(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*tile*/) {return 999;};
    virtual std::tuple<unsigned int, unsigned int, unsigned int> id_to_side_r_phi(unsigned int /*id*/) {return {0, 0, 0};};
    virtual std::tuple<unsigned int, unsigned int, unsigned int> id_to_side_sector_tile(unsigned int /*id*/) {return {0, 0, 0};};
    virtual float r(unsigned int /*id*/) {return -999;};
    virtual float r_from_key(unsigned int /*key*/) {return -999;};
    virtual float r_from_side_r_phi(unsigned int /*side*/, unsigned int /*r_index*/, unsigned int /*phi_index*/) {return -999;};
    virtual float r_from_side_sector_tile(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*tile*/) {return -999;};
    virtual float phi(unsigned int /*id*/) {return -999;};
    virtual float phi_from_key(unsigned int /*key*/) {return -999;};
    virtual float phi_from_side_r_phi(unsigned int /*side*/, unsigned int /*r_index*/, unsigned int /*phi_index*/) {return -999;};
    virtual float phi_from_side_sector_tile(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*tile*/) {return -999;};
    virtual float z(unsigned int /*id*/) {return -999;};
    virtual float z_from_key(unsigned int /*key*/) {return -999;};
    virtual float z_from_side_r_phi(unsigned int /*side*/, unsigned int /*r_inxed*/, unsigned int /*phi_index*/) {return -999;};
    virtual float z_from_side_sector_tile(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*tile*/) {return -999;};
    virtual unsigned int decode_epd(unsigned int /*tower_key*/) {return 999;};

    virtual void Reset() override {return;}; // Reset doesn't need to do anything

  private:
    ClassDefOverride(EpdGeom, 1);
};


#endif // EPD_GEOM_H
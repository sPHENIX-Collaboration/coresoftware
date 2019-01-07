#ifndef PHG4CYLINDERCELL_MVTX_H
#define PHG4CYLINDERCELL_MVTX_H

#include <g4detectors/PHG4CylinderCellv2.h>

#include <cmath>
#include <iostream>
#include <map>

class PHG4CylinderCell_MVTX : public PHG4CylinderCellv2
{
 public:
  PHG4CylinderCell_MVTX();
  virtual ~PHG4CylinderCell_MVTX() {}

  void identify(std::ostream& os = std::cout) const;

  void set_stave_index(const int si) { stave_number = si; }
  int get_stave_index() const { return stave_number; }

  void set_half_stave_index(const int i) { half_stave_number = i; }
  int get_half_stave_index() const { return half_stave_number; }

  void set_module_index(const int i) { module_number = i; }
  int get_module_index() const { return module_number; }

  void set_chip_index(const int i) { chip_number = i; }
  int get_chip_index() const { return chip_number; }

  void set_pixel_index(const int i) { pixel_number = i; }
  int get_pixel_index() const { return pixel_number; }

 protected:
  int stave_number;
  int half_stave_number;
  int module_number;
  int chip_number;
  int pixel_number;

  ClassDef(PHG4CylinderCell_MVTX, 1)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV2_H
#define G4DETECTORS_PHG4CYLINDERCELLV2_H

#include "PHG4CylinderCellv1.h"

#include <iostream>
#include <string>  // for string

class PHG4CylinderCellv2 : public PHG4CylinderCellv1
{
 public:
  PHG4CylinderCellv2();
  ~PHG4CylinderCellv2() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_sensor_index(const std::string& si) override { sensor_index = si; }
  std::string get_sensor_index() const override { return sensor_index; }

  void set_ladder_phi_index(const int i) override { ladder_phi_index = i; }
  int get_ladder_phi_index() const override { return ladder_phi_index; }

  void set_ladder_z_index(const int i) override { ladder_z_index = i; }
  int get_ladder_z_index() const override { return ladder_z_index; }

 protected:
  int ladder_phi_index;
  int ladder_z_index;
  std::string sensor_index;

  ClassDefOverride(PHG4CylinderCellv2, 1)
};

#endif

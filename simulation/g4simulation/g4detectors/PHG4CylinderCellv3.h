// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV3_H
#define G4DETECTORS_PHG4CYLINDERCELLV3_H

#include "PHG4CylinderCellv1.h"

#include <iostream>

class PHG4CylinderCellv3 : public PHG4CylinderCellv1
{
 public:
  PHG4CylinderCellv3();
  ~PHG4CylinderCellv3() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_j_index(const int i) override { j_index = i; }
  int get_j_index() const override { return j_index; }

  void set_k_index(const int i) override { k_index = i; }
  int get_k_index() const override { return k_index; }

  void set_l_index(const int i) override { l_index = i; }
  int get_l_index() const override { return l_index; }

 protected:
  int j_index;
  int k_index;
  int l_index;

  ClassDefOverride(PHG4CylinderCellv3, 1)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALDEFS_H
#define G4DETECTORS_PHG4HCALDEFS_H

#include <string>

namespace PHG4HcalDefs
{
  // parameter names used in lookups by other modules
  static const std::string scipertwr = "n_scinti_plates_per_tower";
  static const std::string innerrad = "inner_radius";
  static const std::string outerrad = "outer_radius";
  static const std::string n_towers = "n_towers";
  static const std::string n_scinti_tiles = "n_scinti_tiles";
  static const std::string n_scinti_tiles_pos = "n_scinti_tiles_pos";
  static const std::string n_scinti_tiles_neg = "n_scinti_tiles_neg";
}  // namespace PHG4HcalDefs

#endif

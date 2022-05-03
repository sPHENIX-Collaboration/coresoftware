// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLDEFS_H
#define G4DETECTORS_PHG4CYLINDERCELLDEFS_H

namespace PHG4CylinderCellDefs
{
  typedef unsigned int keytype;
  static int keybits = 8;
  static int cell_idbits = 32 - keybits;
  enum
  {
    undefined = 0,
    sizebinning = 1,
    etaphibinning = 2,
    etaslatbinning = 3,
    spacalbinning = 4
  };
}  // namespace PHG4CylinderCellDefs

#endif

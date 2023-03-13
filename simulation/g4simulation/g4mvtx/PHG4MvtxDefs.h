// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4MVTX_PHG4MVTXDEFS_H
#define G4MVTX_PHG4MVTXDEFS_H

#include <string>
#include <vector>

namespace PHG4MvtxDefs
{
  static constexpr unsigned int kNLayers = 3;

  enum
  {
    kRmn,
    kRmd,
    kRmx,
    kNModPerStave,
    kPhi0,
    kNStave,
    kNPar
  };

  static const double mvtxdat[kNLayers][kNPar] = {
      {24.61, 25.23, 27.93, 9., 0.285, 12.},  // for each layer: rMin, rMid, rMax, NChip/Stave, phi0, nStaves
      {31.98, 33.36, 36.25, 9., 0.199, 16.},
      {39.93, 41.48, 44.26, 9., 0.166, 20.}};

  static const int GLOBAL = -1;
  static const int ALPIDE_SEGMENTATION = -2;
  static const int SUPPORTPARAMS = -3;

  // passive volume indices

  // detid of support structures

}  // namespace PHG4MvtxDefs

#endif

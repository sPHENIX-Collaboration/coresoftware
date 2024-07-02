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
      {24.610, 25.230, 27.930, 9., 0.2382, 12.},  // for each layer: rMin, rMid, rMax, NChip/Stave, phi0, nStaves
      {31.995, 33.360, 36.258, 9., 0.1937, 16.},
      {39.930, 41.480, 44.260, 9., 0.1481, 20.}};

  static const int GLOBAL = -1;
  static const int ALPIDE_SEGMENTATION = -2;
  static const int SUPPORTPARAMS = -3;

  // passive volume indices

  // detid of support structures

}  // namespace PHG4MvtxDefs

#endif

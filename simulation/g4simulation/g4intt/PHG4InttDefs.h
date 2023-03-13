// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTDEFS_H
#define G4INTT_PHG4INTTDEFS_H

#include <set>

namespace PHG4InttDefs
{
  static const int SEGMENTATION_Z = -1;
  static const int SEGMENTATION_PHI = -2;
  static const int SUPPORTPARAMS = -3;

  static std::set<int> m_SensorSegmentationSet{SEGMENTATION_Z, SEGMENTATION_PHI};

  // passive volume indices
  static const int SI_WAFER = -1;     // strips are distributed on wafer, there are gaps
  static const int SI_INACTIVE = -2;  // inactive si (guard ring etc. I guess)
  static const int HDI_KAPTON = -3;
  static const int HDI_COPPER = -4;
  static const int HDIEXT_KAPTON = -5;
  static const int HDIEXT_COPPER = -6;
  static const int FPHX = -7;
  static const int STAVE_BOTTOM_COOLER = -8;     // changed from PGS
  static const int STAVEEXT_BOTTOM_COOLER = -9;  // changed from PGSEXT
  static const int STAVE_CURVE = -10;
  static const int STAVEEXT_CURVE = -11;
  static const int STAVE_STRAIGHT_OUTER = -12;
  static const int STAVEEXT_STRAIGHT_OUTER = -13;
  static const int STAVE_STRAIGHT_INNER = -14;
  static const int STAVEEXT_STRAIGHT_INNER = -15;
  static const int STAVE_STRAIGHT_COOLER = -16;
  static const int STAVEEXT_STRAIGHT_COOLER = -17;

  static const int SUPPORT_RAIL = -18;
  static const int INTT_OUTER_SKIN = -19;
  static const int INTT_INNER_SKIN = -20;

  static const int SI_GLUE = -21;    // new 29/05/2020
  static const int FPHX_GLUE = -22;  // new 29/05/2020

  static const int SERVICE_BARREL_OUTER = -23;  // new 19/Nov/2021
  static const int SUPPORT_TUBE = -24;          // new 19/Nov/2021
  // detid of support structures
  static const int SUPPORT_DETID = 100;

}  // namespace PHG4InttDefs

#endif

// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SILICONTRACKERDEFS_H
#define G4DETECTORS_PHG4SILICONTRACKERDEFS_H

#include <set>

namespace PHG4SiliconTrackerDefs
{
static const int SEGMENTATION_Z = -1;
static const int SEGMENTATION_PHI = -2;
// this set only exists so we can iterate over the enum
// yes it can be made more fancy but this will do
// and yes stupid CINT does not understand C++11
#ifndef __CINT__
static std::set<int> m_SensorSegmentationSet{SEGMENTATION_Z, SEGMENTATION_PHI};
#endif
// passive volume indices
static const int SI_WAFER = -1;     // strips are distributed on wafer, there are gaps
static const int SI_INACTIVE = -2;  // inactive si (guard ring etc. I guess)
static const int HDI_KAPTON = -3;
static const int HDI_COPPER = -4;
static const int HDIEXT_KAPTON = -5;
static const int HDIEXT_COPPER = -6;
static const int FPHX = -7;
static const int PGS = -8;
static const int PGSEXT = -9;
static const int STAVE_CURVE = -10;
static const int STAVEEXT_CURVE = -11;
static const int STAVE_STRAIGHT_OUTER = -12;
static const int STAVEEXT_STRAIGHT_OUTER = -13;
static const int STAVE_STRAIGHT_INNER = -14;
static const int STAVEEXT_STRAIGHT_INNER = -15;
static const int STAVE_STRAIGHT_COOLER = -16;
static const int STAVEEXT_STRAIGHT_COOLER = -17;
};  // namespace PHG4SiliconTrackerDefs

#endif

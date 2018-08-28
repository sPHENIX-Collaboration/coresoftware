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
};  // namespace PHG4SiliconTrackerDefs

#endif

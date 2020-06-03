// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.

#ifndef MICROMEGASDEFS_H
#define MICROMEGASDEFS_H

/*!
 * \file MicromegasDefs.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <set>

namespace MicromegasDefs
{

  //* tells the direction along which a given cylinder is segmented
  enum class SegmentationType
  {
    SEGMENTATION_Z,
    SEGMENTATION_PHI
  };

}

#endif

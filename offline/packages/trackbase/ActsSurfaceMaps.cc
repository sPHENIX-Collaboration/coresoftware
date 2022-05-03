/*!
 *  \file		ActsSurfaceMaps.cc
 *  \brief		maps hitsetids to Acts Surfaces
 *  \author Tony Frawley <afrawley@fsu.edu>, Joe Osborn <osbornjd@ornl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "ActsSurfaceMaps.h"

#include <Acts/Surfaces/Surface.hpp>

bool ActsSurfaceMaps::isTpcSurface( const Acts::Surface* surface ) const
{ return tpcVolumeIds.find( surface->geometryId().volume() ) != tpcVolumeIds.end(); }
  
bool ActsSurfaceMaps::isMicromegasSurface( const Acts::Surface* surface ) const
{ return micromegasVolumeIds.find( surface->geometryId().volume() ) != micromegasVolumeIds.end(); }

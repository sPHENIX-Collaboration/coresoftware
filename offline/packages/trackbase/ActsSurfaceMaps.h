#ifndef TRACKRECO_ACTSSURFACEMAPS_H
#define TRACKRECO_ACTSSURFACEMAPS_H
/*!
 *  \file		ActsSurfaceMaps.h
 *  \brief		maps hitsetids to Acts Surfaces
 *  \author Tony Frawley <afrawley@fsu.edu>, Joe Osborn <osbornjd@ornl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrkrDefs.h"

namespace Acts{ class Surface; }
class TGeoNode;

#include <map>
#include <memory>
#include <set>
#include <vector>

using Surface = std::shared_ptr<const Acts::Surface>;
using SurfaceVec = std::vector<Surface>;

struct ActsSurfaceMaps
{

  ActsSurfaceMaps() = default;
 
  //! true if given surface corresponds to TPC
  bool isTpcSurface( const Acts::Surface* surface ) const;
    
  //! true if given surface corresponds to Micromegas
  bool isMicromegasSurface( const Acts::Surface* surface ) const;
  
  //! map hitset to Surface for the silicon detectors (MVTX and INTT)
  std::map<TrkrDefs::hitsetkey, Surface> siliconSurfaceMap;

  //! map hitset to surface vector for the TPC
  std::map<unsigned int, SurfaceVec> tpcSurfaceMap;   // uses layer as key

  //! map hitset to surface vector for the micromegas
  std::map<TrkrDefs::hitsetkey, Surface> mmSurfaceMap;
  
  //! map TGeoNode to hitset
  std::map<TrkrDefs::hitsetkey, TGeoNode*> tGeoNodeMap;
 
  //! stores all acts volume ids relevant to the TPC
  /** it is used to quickly tell if a given Acts Surface belongs to the TPC */
  std::set<int> tpcVolumeIds;

  //! stores all acts volume ids relevant to the micromegas
  /** it is used to quickly tell if a given Acts Surface belongs to micromegas */
  std::set<int> micromegasVolumeIds;

};

#endif

#ifndef TRACKBASE_ACTSSURFACEMAPS_H
#define TRACKBASE_ACTSSURFACEMAPS_H
/*!
 *  \file		ActsSurfaceMaps.h
 *  \brief		maps hitsetids to Acts Surfaces
 *  \author Tony Frawley <afrawley@fsu.edu>, Joe Osborn <osbornjd@ornl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrkrDefs.h"
#include "ActsTrackingGeometry.h"

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

namespace Acts{ class Surface; }
class TGeoNode;
class TrkrCluster;

#include <map>
#include <memory>
#include <set>
#include <vector>

using Surface = std::shared_ptr<const Acts::Surface>;
using SurfaceVec = std::vector<Surface>;

struct ActsSurfaceMaps
{
 public:
  ActsSurfaceMaps() = default;
 
  //! true if given surface corresponds to TPC
  bool isTpcSurface( const Acts::Surface* surface ) const;
    
  //! true if given surface corresponds to Micromegas
  bool isMicromegasSurface( const Acts::Surface* surface ) const;
  
  Surface getSurface(TrkrDefs::cluskey, TrkrCluster* cluster) const;

  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey) const;

  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey,
    TrkrDefs::subsurfkey surfkey) const;

  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey) const;

  //! map hitset to Surface for the silicon detectors (MVTX and INTT)
  std::map<TrkrDefs::hitsetkey, Surface> m_siliconSurfaceMap;

  //! map hitset to surface vector for the TPC
  std::map<unsigned int, SurfaceVec> m_tpcSurfaceMap;   // uses layer as key

  //! map hitset to surface vector for the micromegas
  std::map<TrkrDefs::hitsetkey, Surface> m_mmSurfaceMap;
  
  //! map TGeoNode to hitset
  std::map<TrkrDefs::hitsetkey, TGeoNode*> m_tGeoNodeMap;
 
 //! stores all acts volume ids relevant to the TPC
  /** it is used to quickly tell if a given Acts Surface belongs to the TPC */
  std::set<int> m_tpcVolumeIds;

  //! stores all acts volume ids relevant to the micromegas
  /** it is used to quickly tell if a given Acts Surface belongs to micromegas */
  std::set<int> m_micromegasVolumeIds;

};

#endif

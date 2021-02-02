#ifndef TRACKRECO_ACTSSURFACEMAPS_H
#define TRACKRECO_ACTSSURFACEMAPS_H

#include <trackbase/TrkrDefs.h>
#include <TGeoNode.h>

#include <map>
#include <vector>

using Surface = std::shared_ptr<const Acts::Surface>;
using SurfaceVec = std::vector<Surface>;


struct ActsSurfaceMaps
{

  ActsSurfaceMaps(){};
  std::map<TrkrDefs::hitsetkey, Surface> siliconSurfaceMap;
  std::map<TrkrDefs::hitsetkey, SurfaceVec> tpcSurfaceMap;
  std::map<TrkrDefs::hitsetkey, SurfaceVec> mmSurfaceMap;
  std::map<TrkrDefs::hitsetkey, TGeoNode*> tGeoNodeMap;
 
};

#endif

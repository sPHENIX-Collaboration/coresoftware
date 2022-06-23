#include "ActsGeometry.h"
#include "TrkrCluster.h"

namespace
{
  /// square
  template<class T> inline constexpr T square(const T& x) {return x*x;}

  /// get radius from coordinates
  template<class T> T radius(const T& x, const T& y)
  { return std::sqrt(square(x) + square(y));}
}

Eigen::Matrix<float,3,1> ActsGeometry::getGlobalPositionF(
  TrkrDefs:: cluskey key,       
  TrkrCluster* cluster) const
{
  const Acts::Vector3 doublePos = getGlobalPosition(key, cluster);
  return Eigen::Matrix<float,3,1>(doublePos(0), doublePos(1), doublePos(2));
}

Acts::Vector3 ActsGeometry::getGlobalPosition(TrkrDefs:: cluskey key,
  TrkrCluster* cluster) const
{
  Acts::Vector3 glob;
 
  const auto trkrid = TrkrDefs::getTrkrId(key);

  auto surface = maps().getSurface(key, cluster);

  if(!surface)
    {
   
      std::cerr << "Couldn't identify cluster surface. Returning NAN"
		<< std::endl;
      glob(0) = NAN;
      glob(1) = NAN;
      glob(2) = NAN;
      return glob;
    } 

  Acts::Vector2 local(cluster->getLocalX(), cluster->getLocalY());
  Acts::Vector3 global;
  /// If silicon/TPOT, the transform is one-to-one since the surface is planar
  if(trkrid != TrkrDefs::tpcId)
    {
      global = surface->localToGlobal(geometry().geoContext,
				      local * Acts::UnitConstants::cm,
				      Acts::Vector3(1,1,1));
      global /= Acts::UnitConstants::cm;
      return global;
    }

  /// Otherwise do the manual calculation
  /// Undo the manual calculation that is performed in TpcClusterizer
  auto surfCenter = surface->center(geometry().geoContext);

  surfCenter /= Acts::UnitConstants::cm;
  double surfPhiCenter = atan2(surfCenter(1), surfCenter(0));
  double surfRadius = radius(surfCenter(0), surfCenter(1));
  double surfRPhiCenter = surfPhiCenter * surfRadius;
  
  double clusRPhi = local(0) + surfRPhiCenter;
  double gclusz = local(1) + surfCenter(2);
  
  double clusphi = clusRPhi / surfRadius;
  double gclusx = surfRadius * cos(clusphi);
  double gclusy = surfRadius * sin(clusphi);

  global(0) = gclusx;
  global(1) = gclusy;
  global(2) = gclusz;
  
  return global;
}


Surface ActsGeometry::get_tpc_surface_from_coords(
  TrkrDefs::hitsetkey hitsetkey,
  Acts::Vector3 world,
  TrkrDefs::subsurfkey& subsurfkey) const
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  auto tpcmap = maps().m_tpcSurfaceMap;
  auto mapIter = tpcmap.find(layer);
  
  if(mapIter == tpcmap.end())
    {
      std::cout << "Error: hitsetkey not found in ActsGeometry::get_tpc_surface_from_coords, hitsetkey = "
		<< hitsetkey << std::endl;
      return nullptr;
    }
  
  double world_phi = atan2(world[1], world[0]);
  double world_z = world[2];
  
  std::vector<Surface> surf_vec = mapIter->second;
  unsigned int surf_index = 999;
    
  // Predict which surface index this phi and z will correspond to
  // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
  double fraction =  (world_phi + M_PI) / (2.0 * M_PI);
  double rounded_nsurf = round( (double) (surf_vec.size()/2) * fraction  - 0.5);
  unsigned int nsurf = (unsigned int) rounded_nsurf; 
  if(world_z < 0)
    { nsurf += surf_vec.size()/2; }

  Surface this_surf = surf_vec[nsurf];
      
  auto vec3d = this_surf->center(geometry().geoContext);
  std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
  double surf_z = surf_center[2];
  double surf_phi = atan2(surf_center[1], surf_center[0]);
  double surfStepPhi = geometry().tpcSurfStepPhi;
  double surfStepZ = geometry().tpcSurfStepZ;

  if( (world_phi > surf_phi - surfStepPhi / 2.0 && world_phi < surf_phi + surfStepPhi / 2.0 ) &&
      (world_z > surf_z - surfStepZ / 2.0 && world_z < surf_z + surfStepZ / 2.0) )	
    {
      surf_index = nsurf;
      subsurfkey = nsurf;
    }    
  else
    {
      return nullptr;
    }
  
  return surf_vec[surf_index];
  
}


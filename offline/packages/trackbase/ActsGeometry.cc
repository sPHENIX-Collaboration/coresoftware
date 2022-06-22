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

/*!
 *  \file		ActsSurfaceMaps.cc
 *  \brief		maps hitsetids to Acts Surfaces
 *  \author Tony Frawley <afrawley@fsu.edu>, Joe Osborn <osbornjd@ornl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "ActsSurfaceMaps.h"
#include "TrkrCluster.h"
#include "MvtxDefs.h"
#include "InttDefs.h"

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

namespace
{
  /// square
  template<class T> inline constexpr T square(const T& x) {return x*x;}

  /// get radius from coordinates
  template<class T> T radius(const T& x, const T& y)
  { return std::sqrt(square(x) + square(y));}
}


bool ActsSurfaceMaps::isTpcSurface( const Acts::Surface* surface ) const
{ return m_tpcVolumeIds.find( surface->geometryId().volume() ) != m_tpcVolumeIds.end(); }
  
bool ActsSurfaceMaps::isMicromegasSurface( const Acts::Surface* surface ) const
{ return m_micromegasVolumeIds.find( surface->geometryId().volume() ) != m_micromegasVolumeIds.end(); }



Surface ActsSurfaceMaps::getSurface(TrkrDefs::cluskey key,       
					TrkrCluster *cluster) const
{
  const auto trkrid = TrkrDefs::getTrkrId(key);
  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  switch( trkrid )
  {
    case TrkrDefs::TrkrId::micromegasId: 
    {
      return getMMSurface( hitsetkey);
    }
    
    case TrkrDefs::TrkrId::tpcId: 
    {
      const auto surfkey = cluster->getSubSurfKey();
      return getTpcSurface(hitsetkey, surfkey);
    }
    
    case TrkrDefs::TrkrId::mvtxId:
    case TrkrDefs::TrkrId::inttId:
    {
      return getSiliconSurface(hitsetkey);
    }
  }
  
  // unreachable
  return nullptr;
}

Surface ActsSurfaceMaps::getSiliconSurface(TrkrDefs::hitsetkey hitsetkey) const
{
  unsigned int trkrid =  TrkrDefs::getTrkrId(hitsetkey);
  TrkrDefs::hitsetkey tmpkey = hitsetkey;

  if(trkrid == TrkrDefs::inttId)
    {
      // Set the hitsetkey crossing to zero
      tmpkey = InttDefs::resetCrossingHitSetKey(hitsetkey);
    }

  if(trkrid == TrkrDefs::mvtxId)
    {
      // Set the hitsetkey crossing to zero
      tmpkey = MvtxDefs::resetStrobeHitSetKey(hitsetkey);
    }
     
  auto iter = m_siliconSurfaceMap.find(tmpkey);
  if(iter != m_siliconSurfaceMap.end())
    {
      return iter->second;
    }
  
  /// If it can't be found, return nullptr
  {
    std::cout << "Failed to find silicon surface for hitsetkey " << hitsetkey << " tmpkey " << tmpkey << std::endl;
  }  
  return nullptr;
}

Surface ActsSurfaceMaps::getTpcSurface(TrkrDefs::hitsetkey hitsetkey, 
				       TrkrDefs::subsurfkey surfkey) const
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const auto iter = m_tpcSurfaceMap.find(layer);
  
  if(iter != m_tpcSurfaceMap.end())
  {
    auto surfvec = iter->second;
    return surfvec.at(surfkey);
  }

  /// If it can't be found, return nullptr to skip this cluster
  return nullptr;
}


Surface ActsSurfaceMaps::getMMSurface(TrkrDefs::hitsetkey hitsetkey) const
{
  const auto iter = m_mmSurfaceMap.find( hitsetkey );
  return (iter == m_mmSurfaceMap.end()) ? nullptr:iter->second;
}


Eigen::Matrix<float,3,1> ActsSurfaceMaps::getGlobalPositionF(
  TrkrDefs:: cluskey key,       
  TrkrCluster* cluster,
  ActsTrackingGeometry* tGeometry) const
{
  const Acts::Vector3 doublePos = getGlobalPosition(key, cluster, tGeometry);
  return Eigen::Matrix<float,3,1>(doublePos(0), doublePos(1), doublePos(2));
}

Acts::Vector3 ActsSurfaceMaps::getGlobalPosition(
  TrkrDefs:: cluskey key,
  TrkrCluster* cluster,
  ActsTrackingGeometry *tGeometry) const
{

  Acts::Vector3 glob;
 
  const auto trkrid = TrkrDefs::getTrkrId(key);

  auto surface = getSurface(key, cluster);

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
      global = surface->localToGlobal(tGeometry->geoContext,
				      local * Acts::UnitConstants::cm,
				      Acts::Vector3(1,1,1));
      global /= Acts::UnitConstants::cm;
      return global;
    }

  /// Otherwise do the manual calculation
  /// Undo the manual calculation that is performed in TpcClusterizer
  auto surfCenter = surface->center(tGeometry->geoContext);

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

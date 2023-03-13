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

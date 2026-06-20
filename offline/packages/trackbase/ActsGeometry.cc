#include "ActsGeometry.h"
#include "TpcDefs.h"
#include "TrkrCluster.h"
#include "alignmentTransformationContainer.h"

#include <phool/sphenix_constants.h>

#include <Acts/Definitions/Algebra.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

namespace
{
  /// square
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  /// get radius from coordinates
  template <class T>
  T radius(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }
}  // namespace

//________________________________________________________________________________________________
Acts::Vector3 ActsGeometry::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  Acts::Vector3 glob;

  const auto trkrid = TrkrDefs::getTrkrId(key);
  if (trkrid == TrkrDefs::tpcId)
  {
    return getGlobalPositionTpc(key, cluster);
  }

  /// If silicon/TPOT, the transform is one-to-one since the surface is planar

  auto surface = m_surfMaps.getSurface(key, cluster);

  if (!surface)
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
  global = surface->localToGlobal(m_tGeometry.getGeoContext(),
                                  local * Acts::UnitConstants::cm,
                                  Acts::Vector3(1, 1, 1));
  global /= Acts::UnitConstants::cm;

  return global;
}

//________________________________________________________________________________________________
Acts::Vector3 ActsGeometry::getGlobalPositionTpc(const TrkrDefs::hitsetkey& hitsetkey,
const TrkrDefs::hitkey& hitkey, const float& phi, const float& rad,
const float& clockPeriod) const
{
  Acts::Vector3 glob;
  const auto trkrid = TrkrDefs::getTrkrId(hitsetkey);
  if (trkrid != TrkrDefs::tpcId)
  {
    std::cout << "ActsGeometry::getGlobalPositionTpc -  this is the wrong global transform for silicon or MM's clusters! Returning zero" << std::endl;
    return glob;
  }

  auto tbin = TpcDefs::getTBin(hitkey);
  double zdriftlength = tbin * clockPeriod * _drift_velocity;  // cm
  double zloc =  _max_driftlength / 2.0 - zdriftlength;                   // local z relative to surface center (for north side):
  unsigned int side = TpcDefs::getSide(hitsetkey);
  if (side == 0)
  {
    zloc = -zloc;
  }
  
  float surfaceZCenter = _max_driftlength/2.0 + _CM_halfwidth;
  
  auto x = rad * std::cos(phi);
  auto y = rad * std::sin(phi);
  auto z = surfaceZCenter + zloc;
  glob.x() = x;
  glob.y() = y;
  glob.z() = z;
  return glob;
}

//________________________________________________________________________________________________
Acts::Vector3 ActsGeometry::getGlobalPositionTpc(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  Acts::Vector3 glob;

  // This method is for the TPC only
  const auto trkrid = TrkrDefs::getTrkrId(key);
  if (trkrid != TrkrDefs::tpcId)
  {
    std::cout << "ActsGeometry::getGlobalPositionTpc -  this is the wrong global transform for silicon or MM's clusters! Returning zero" << std::endl;
    return glob;
  }

  auto surface = m_surfMaps.getSurface(key, cluster);

  /*
  std::cout << " getGlobalPositionTpc transform is: " << std::endl
	    <<  surface->transform(m_tGeometry.getGeoContext()).matrix()
	    << std::endl;
  alignmentTransformationContainer::use_alignment = false;
  Acts::Vector3 ideal_center = surface->center(m_tGeometry.getGeoContext()) * 0.1;
  alignmentTransformationContainer::use_alignment = true;  
  Acts::Vector3 sensorCenter = surface->center(m_tGeometry.getGeoContext()) * 0.1;  // cm
  std::cout << "  ideal surface center: " << ideal_center << std::endl;
  std::cout << "  aligned surface center: " << sensorCenter << std::endl;
  */
  
  if (!surface)
  {
    std::cerr << "Couldn't identify cluster surface. Returning NAN"
              << std::endl;
    glob(0) = NAN;
    glob(1) = NAN;
    glob(2) = NAN;
    return glob;
  }

  Acts::Vector2 local = getLocalCoords(key, cluster);  // no crossing correction here
  
  glob = surface->localToGlobal(m_tGeometry.getGeoContext(),
                                local * Acts::UnitConstants::cm,
                                Acts::Vector3(1, 1, 1));
  glob /= Acts::UnitConstants::cm;


  // std::cout << "  local " << local << std::endl;
  // std::cout << "  glob " << glob << std::endl;

  
  return glob;
}

Surface ActsGeometry::get_tpc_surface_from_coords(
    TrkrDefs::hitsetkey hitsetkey,
    Acts::Vector3 world,
    TrkrDefs::subsurfkey& subsurfkey) const
{
  // Assume that the world coordinates are in the sPHENIX frame, where the TPC is tilted
  // We convert the position to tpc envelope coordinates, where we know where everything is
  Acts::Vector3 world_envelope = transformTpcWorldToEnvelope(world);
  double world_phi = atan2(world_envelope[1], world_envelope[0]);

  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int side = TpcDefs::getSide(hitsetkey);
  unsigned int sector = TpcDefs::getSectorId(hitsetkey);
 
  // returns an iterator to all of the surfaces for this layer
  auto mapIter = m_surfMaps.m_tpcSurfaceMap.find(layer);

  if (mapIter == m_surfMaps.m_tpcSurfaceMap.end())
  {
    std::cout << "Error: hitsetkey not found in ActsGeometry::get_tpc_surface_from_coords, hitsetkey = "
              << hitsetkey << std::endl;
    return nullptr;
  }

  const auto& surf_vec = mapIter->second;
  unsigned int surf_index = 999;

  // Apparently, tilting the TPC leads to the surfaces not being sorted in phi in the outer layers
  // just test all surfaces in each layer for now
  for(unsigned int isurf = 0; isurf < surf_vec.size(); ++isurf)
    {
      Surface this_surf = surf_vec[isurf];
      auto surf_center = this_surf->center(m_tGeometry.getGeoContext());
      surf_center /= 10.0;   // convert from mm to cm
      //this surface center includes the TPC tilt used in PHG4TpcDetector construction, transform it to tpc envelope coordinates
      Acts::Vector3 surf_center_envelope = transformTpcWorldToEnvelope(surf_center);
      double surf_phi = atan2(surf_center_envelope[1], surf_center_envelope[0]);  
      double surfStepPhi = m_tGeometry.tpcSurfStepPhi;

      if ((world_phi > surf_phi - surfStepPhi / 2.0) && (world_phi < surf_phi + surfStepPhi / 2.0))
	{
	  if(surf_center.z() < 0 && side != 0) { continue; }
	  if(surf_center.z() > 0 && side != 1) { continue; }
	  surf_index = isurf;
	  subsurfkey = isurf;
	  break;
	}            
    }  

  if(surf_index == 999)
    {
    std::cout << "Error: surface not found in ActsGeometry::get_tpc_surface_from_coords "
              << " layer " << layer << " side " << side << " sector " << sector << " world_phi " << world_phi << "  world[0]  " << world[0] << " world[1] " << world[1] << " hitsetkey " << hitsetkey << std::endl;
    return nullptr;
    }

  return surf_vec[surf_index];
  
  /*
  // Predict which surface index this phi and side will correspond to
  // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
  // we use TPC side from the hitsetkey, since z can be either sign in north and south, depending on crossing
  double fraction = (world_phi + M_PI) / (2.0 * M_PI);

  double rounded_nsurf = std::round((double) (surf_vec.size() / 2) * fraction - 0.5);  // NOLINT
  unsigned int nsurfm = (unsigned int) rounded_nsurf;
 std::cout << " surf_vec.size " << surf_vec.size() << " rounded_nsurf " << rounded_nsurf << " initial nsurfm " << nsurfm << std::endl;
  
  if (side == 0)
  {
    nsurfm += surf_vec.size() / 2;
  }
  unsigned int nsurf = nsurfm % surf_vec.size();
  Surface this_surf = surf_vec[nsurf];
  auto surf_center = this_surf->center(m_tGeometry.getGeoContext());
  surf_center /= 10.0;   // convert from mm to cm
  //this surface center is from the default geometry, which includes the TPC tilt used in PHG4TpcDetector construction
  // transform it to tpc envelope coordinates
  Acts::Vector3 surf_center_envelope = m_tpc_world_envelope_transform * surf_center;

  double surf_phi = atan2(surf_center_envelope[1], surf_center_envelope[0]);  
  double surfStepPhi = m_tGeometry.tpcSurfStepPhi;
  
  if ((world_phi > surf_phi - surfStepPhi / 2.0) && (world_phi < surf_phi + surfStepPhi / 2.0))
  {
    surf_index = nsurf;
    subsurfkey = nsurf;
    std::cout << "success, found nsurf = " << nsurf << std::endl;
  }
  else
  {
    // check for the periodic boundary condition
    auto firstsurf = *surf_vec.begin();
    auto firstsurfcenter = firstsurf->center(m_tGeometry.getGeoContext());
    firstsurfcenter /= 10.0;
    auto firstsurfcenter_envelope = m_tpc_world_envelope_transform * firstsurfcenter;
    double firstsurf_phi = atan2(firstsurfcenter_envelope[1], firstsurfcenter_envelope[0]);
    if (world_phi < -M_PI)
    {
      world_phi += 2.0 * M_PI;
    }
    //check a few surfaces around this one
    for( int i = -1; i <= 1; i++)
    {
      if(i==0) // already tried this one
      {
        continue;
      }
      unsigned int new_nsurf = (nsurf+i) % surf_vec.size();
      this_surf = surf_vec[new_nsurf];
      surf_center = this_surf->center(m_tGeometry.getGeoContext());
      surf_center /= 10.0;
      surf_center_envelope = m_tpc_world_envelope_transform * surf_center;
      surf_phi = atan2(surf_center_envelope[1], surf_center_envelope[0]);
      double this_philow =  surf_phi - surfStepPhi / 2.0;
      double this_phihigh =  surf_phi + surfStepPhi / 2.0;
      
      if ((world_phi > this_philow) && (world_phi < this_phihigh))
      {
	std::cout << "success, found nsurf = " << new_nsurf << std::endl;	
        surf_index = new_nsurf;
        subsurfkey = new_nsurf;
        return surf_vec[surf_index];
      }
    }
    return nullptr;
  }
  */
  

}

//________________________________________________________________________________________________
Acts::Transform3 ActsGeometry::makeAffineTransform(Acts::Vector3 rot, Acts::Vector3 trans) const
{
  Acts::Transform3 actsAffine;

  Eigen::AngleAxisd alpha(rot(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd beta(rot(1), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd gamma(rot(2), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> q = gamma * beta * alpha;
  actsAffine.linear() = q.matrix();

  Eigen::Vector3d translation(trans(0), trans(1), trans(2));
  actsAffine.translation() = translation;

  return actsAffine;
}

//________________________________________________________________________________________________
Acts::Vector2 ActsGeometry::getLocalCoords(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  short int crossing = 0;
  Acts::Vector2 local = getLocalCoords(key, cluster, crossing);
  return local;
}

  
//________________________________________________________________________________________________
Acts::Vector2 ActsGeometry::getLocalCoords(TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing) const
{
  Acts::Vector2 local;

  const auto trkrid = TrkrDefs::getTrkrId(key);
  if (trkrid == TrkrDefs::tpcId)
  {
    double crossing_tzero_correction = crossing * sphenix_constants::time_between_crossings;
    double tcorrected = cluster->getLocalY() +  _tpc_tzero + _sampa_tzero_bias - crossing_tzero_correction;
    double zdriftlength = tcorrected * _drift_velocity; 
    double zloc = _max_driftlength/2.0 - zdriftlength;         // local z relative to surface center (for north side):
    unsigned int side = TpcDefs::getSide(key);
    if (side == 0)
    {
      zloc = -zloc;
    }
    local(0) = cluster->getLocalX();
    local(1) = zloc;

    /*
    std::cout << " clust " << cluster->getLocalY() << " tpc tzero " << _tpc_tzero << " sampa tbias " << _sampa_tzero_bias
	      << " crossing tzero correction " << crossing_tzero_correction << " corrected clust " << tcorrected
	      << " drift vel " << _drift_velocity 
	      << " crossing " << crossing << " crossing period " << sphenix_constants::time_between_crossings
	      << " maxdriftlength " << _max_driftlength << " zdriftlength " << zdriftlength << " zloc " << zloc << std::endl;
    */
    
  }
  else
  {
    local(0) = cluster->getLocalX();
    local(1) = cluster->getLocalY();
  }

  return local;
}

  Acts::Vector3  ActsGeometry::transformTpcWorldToEnvelope(Acts::Vector3 world) const
  {
    Acts::Vector3 envelope = m_tpc_world_envelope_transform * world;

    return envelope;
  }

  Acts::Vector3  ActsGeometry::transformTpcEnvelopeToWorld(Acts::Vector3 envelope) const
  {
    Acts::Vector3 world = m_tpc_world_envelope_transform.inverse() * envelope;

    return world;
  }

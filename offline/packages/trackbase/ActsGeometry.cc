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

  /*
    // Leaving this here for now, since it gives comprehensive diagnostic info 

    unsigned int sskey = cluster->getSubSurfKey();
    unsigned int layer = TrkrDefs::getLayer(key);
    unsigned int side = TpcDefs::getSide(key);
    unsigned int sector = TpcDefs::getSectorId(key);

    alignmentTransformationContainer::use_alignment = false;
    Acts::Vector3 ideal_center = surface->center(m_tGeometry.getGeoContext()) * 0.1;
    alignmentTransformationContainer::use_alignment = true;  
    
    Acts::Vector3 env_center = m_tpc_world_envelope_transform * ideal_center;
    Acts::Vector3 sensorCenter = surface->center(m_tGeometry.getGeoContext()) * 0.1;  // cm

    std::cout << "  layer " << layer << " side " << side << " sector " << sector << " sskey " << sskey
	    << "  aligned global: " << glob[0] << "  " << glob[1] << "  " << glob[2]  
	    << "  aligned surf: " << sensorCenter[0] << "  " << sensorCenter[1] << "  " << sensorCenter[2]
    	    << "  ideal surf: " << ideal_center[0] << "  " << ideal_center[1] << "  " << ideal_center[2]
	    << "  env surf: " << env_center[0] << "  " << env_center[1] << "  " << env_center[2]
	    << std::endl;
  */
  
  return glob;
}

Surface ActsGeometry::get_tpc_surface_from_coords(
    TrkrDefs::hitsetkey hitsetkey,
    Acts::Vector3 cluster,
    TrkrDefs::subsurfkey& subsurfkey) const
{
  // This method finds the subsurface for global cluster positions

  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int side = TpcDefs::getSide(hitsetkey);
  unsigned int sector = TpcDefs::getSectorId(hitsetkey);

  //    std::cout << " get_tpc_surface_from_coords: layer " << layer << " side " << side << " sector " << sector << std::endl;
    
  double surfStepPhi = m_tGeometry.tpcSurfStepPhi;
	
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

  // convert global to mm for consistency with surface centers
  cluster *= 10.0;
  
  // Apparently, tilting the TPC leads to the surfaces not being sorted in phi in some layers
  // just test all surfaces in each layer for now
  double min_dphi = 999.0;
  unsigned int min_surf_index = 999;
  for(unsigned int isurf = 0; isurf < surf_vec.size(); ++isurf)
    {
      Surface this_surf = surf_vec[isurf];

      // quick check to eliminate wrong side surfaces 
      auto surf_center_test = this_surf->center(m_tGeometry.getGeoContext());
      unsigned int surf_side = 0;
      if(surf_center_test.z() > 0.0) { surf_side = 1; }
      if(surf_side != side) { continue; }
      
      // the cluster coordinates are in the sPHENIX frame, where the TPC is tilted and alignment
      // transforms are implemented. We must convert the cluster  to tpc envelope coordinates,
      // where we know where the fake surfaces are located.
      //    so, transform:  cluster_aligned->local->cluster_noalign->envelope
      
      Acts::Vector3 local = this_surf->localToGlobalTransform(m_tGeometry.getGeoContext()).inverse() * (cluster);
      // transform local to unaligned geometry (simulation) global
      alignmentTransformationContainer::use_alignment = false;
      Acts::Vector3 cluster_noalign = this_surf->localToGlobalTransform(m_tGeometry.getGeoContext()) * (local);   
      auto surf_center_noalign = this_surf->center(m_tGeometry.getGeoContext());
      alignmentTransformationContainer::use_alignment = true;

      // transform simulation geometry global to envelope coords
      Acts::Vector3 cluster_envelope = transformTpcWorldToEnvelope(cluster_noalign/10.0) * 10.0;  // transform needs cm
      Acts::Vector3 surf_center_envelope = transformTpcWorldToEnvelope(surf_center_noalign/10.0) * 10.0;  
      
      double cluster_phi_envelope = atan2(cluster_envelope[1], cluster_envelope[0]);
      double surf_phi_envelope = atan2(surf_center_envelope[1], surf_center_envelope[0]);
      const double dphi = std::atan2(std::sin(cluster_phi_envelope - surf_phi_envelope), std::cos(cluster_phi_envelope - surf_phi_envelope));
      
      if(std::abs(dphi) < min_dphi)
	{
	  min_dphi = std::abs(dphi);
	  min_surf_index = isurf;
	}
    }
  
  surf_index = min_surf_index;
  subsurfkey = min_surf_index;
  
  if(min_dphi > surfStepPhi)
    {
      // too large to be due to cluster uncertainty
      std::cout << "Error: surface not found in ActsGeometry::get_tpc_surface_from_coords "
		<< " layer " << layer << " side " << side << " sector " << sector
		<< " min_dphi " << min_dphi << " min_surf_index " << min_surf_index
		<< " cluster[0]  " << cluster[0] << " cluster[1] " << cluster[1] << " cluster[2] " << cluster[2] << " mm "
		<< " cluster phi " << atan2(cluster[1],cluster[0])
		<< " hitsetkey " << hitsetkey << std::endl;
      return nullptr;
    }
  
  return surf_vec[surf_index];
}

Surface ActsGeometry::get_clusterizer_tpc_surface(
    TrkrDefs::hitsetkey hitsetkey,
    Acts::Vector3 clus_envelope,
    TrkrDefs::subsurfkey& subsurfkey) const
{
  // This method finds the subsurface for the clusterizer so the key can be added to the cluster
  // The input cluster position is in envelope coordinates

  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int side = TpcDefs::getSide(hitsetkey);
  unsigned int sector = TpcDefs::getSectorId(hitsetkey);

  double surfStepPhi = m_tGeometry.tpcSurfStepPhi;
  
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

  // convert position to mm for consistency with surfaces
  clus_envelope *= 10.0;

  double clus_phi_envelope = atan2(clus_envelope[1], clus_envelope[0]);
      
  // Apparently, tilting the TPC leads to the surfaces not being sorted in phi in the outer layers
  // just test all surfaces in each layer for now
  double min_dphi = 999.0;
  unsigned int min_surf_index = 999;
  for(unsigned int isurf = 0; isurf < surf_vec.size(); ++isurf)
    {
      Surface this_surf = surf_vec[isurf];
      
      // get the surface center before alignment
      // leave the alignment flag the way you found it!
      Acts::Vector3 surf_center_noalign(0,0,0);
      Acts::Vector3 surf_center_local(0,0,0);
      if(alignmentTransformationContainer::use_alignment)
	{
	  alignmentTransformationContainer::use_alignment = false;
	  surf_center_noalign = this_surf->center(m_tGeometry.getGeoContext());
	  surf_center_local = this_surf->localToGlobalTransform(m_tGeometry.getGeoContext()).inverse() * (surf_center_noalign);
	  alignmentTransformationContainer::use_alignment = true;	  
	}
      else
	{
	  // this will be the case when called from the clusterizer
	  surf_center_noalign = this_surf->center(m_tGeometry.getGeoContext());
	  surf_center_local = this_surf->localToGlobalTransform(m_tGeometry.getGeoContext()).inverse() * (surf_center_noalign);
	}

      // eliminate wrong side surfaces 
      unsigned int surf_side = 0;
      if(surf_center_noalign.z() > 0.0) { surf_side = 1; }
      if(surf_side != side) { continue; }

      
      // transform simulation geometry global to envelope coords
      Acts::Vector3 surf_center_envelope = transformTpcWorldToEnvelope(surf_center_noalign/10.0) * 10.0;    // transform uses cm
      double surf_phi_envelope = atan2(surf_center_envelope[1], surf_center_envelope[0]);  
      const double dphi = std::atan2(std::sin(clus_phi_envelope - surf_phi_envelope), std::cos(clus_phi_envelope - surf_phi_envelope));

      if(std::abs(dphi) < min_dphi)
	{
	  min_dphi = std::abs(dphi);
	  min_surf_index = isurf;
	}
    }
  
  surf_index = min_surf_index;
  subsurfkey = min_surf_index;
  
  if(min_dphi > surfStepPhi)
    {
      // too large to be due to cluster uncertainty
      std::cout << "Error: surface not found in ActsGeometry::get_tpc_surface_from_coords "
		<< " layer " << layer << " side " << side << " sector " << sector
		<< " min_dphi " << min_dphi << " min_surf_index " << min_surf_index
		<< " cluster[0]  " << clus_envelope[0] << " cluster[1] " << clus_envelope[1] << " cluster[2] " << clus_envelope[2] << " mm "
		<< " cluster phi " << atan2(clus_envelope[1],clus_envelope[0])
		<< " hitsetkey " << hitsetkey << std::endl;
      return nullptr;
    }
  
  return surf_vec[surf_index];
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

  Acts::Vector3  ActsGeometry::transformTpcWorldToEnvelope(const Acts::Vector3& world) const
  {
    Acts::Vector3 envelope = m_tpc_world_envelope_transform * world;

    return envelope;
  }

  Acts::Vector3  ActsGeometry::transformTpcEnvelopeToWorld(const Acts::Vector3& envelope) const
  {
    Acts::Vector3 world = m_tpc_world_envelope_transform.inverse() * envelope;

    return world;
  }

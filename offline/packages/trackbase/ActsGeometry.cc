#include "ActsGeometry.h"
#include <Acts/Definitions/Algebra.hpp>
#include "TpcDefs.h"
#include "TrkrCluster.h"
#include "alignmentTransformationContainer.h"
#include <phool/sphenix_constants.h>

namespace
{
  /// square
  template <class T>
  inline constexpr T square(const T& x)
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

  return glob;
}

Surface ActsGeometry::get_tpc_surface_from_coords(
    TrkrDefs::hitsetkey hitsetkey,
    Acts::Vector3 world,
    TrkrDefs::subsurfkey& subsurfkey) const
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int side = TpcDefs::getSide(hitsetkey);

  auto mapIter = m_surfMaps.m_tpcSurfaceMap.find(layer);

  if (mapIter == m_surfMaps.m_tpcSurfaceMap.end())
  {
    std::cout << "Error: hitsetkey not found in ActsGeometry::get_tpc_surface_from_coords, hitsetkey = "
              << hitsetkey << std::endl;
    return nullptr;
  }
  double world_phi = atan2(world[1], world[0]);

  const auto& surf_vec = mapIter->second;
  unsigned int surf_index = 999;

  // Predict which surface index this phi and side will correspond to
  // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
  // we use TPC side from the hitsetkey, since z can be either sign in northa nd south, depending on crossing
  double fraction = (world_phi + M_PI) / (2.0 * M_PI);

  double rounded_nsurf = std::round((double) (surf_vec.size() / 2) * fraction - 0.5);  // NOLINT
  unsigned int nsurfm = (unsigned int) rounded_nsurf;

  if (side == 0)
  {
    nsurfm += surf_vec.size() / 2;
  }

  unsigned int nsurf = nsurfm % surf_vec.size();
  Surface this_surf = surf_vec[nsurf];

  auto vec3d = this_surf->center(m_tGeometry.getGeoContext());
  std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
  double surf_phi = atan2(surf_center[1], surf_center[0]);
  double surfStepPhi = m_tGeometry.tpcSurfStepPhi;

  if ((world_phi > surf_phi - surfStepPhi / 2.0 && world_phi < surf_phi + surfStepPhi / 2.0))
  {
    surf_index = nsurf;
    subsurfkey = nsurf;
  }
  else
  {
    // check for the periodic boundary condition
    auto firstsurf = *surf_vec.begin();
    auto firstsurfcenter = firstsurf->center(geometry().getGeoContext());
    float firstsurf_phi = atan2(firstsurfcenter[1], firstsurfcenter[0]);
    if (world_phi < firstsurf_phi - surfStepPhi / 2.0)
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
      vec3d = this_surf->center(geometry().getGeoContext());
      surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
      surf_phi = atan2(surf_center[1], surf_center[0]);

      if ((world_phi > surf_phi - surfStepPhi / 2.0 && world_phi < surf_phi + surfStepPhi / 2.0))
      {
        surf_index = new_nsurf;
        subsurfkey = new_nsurf;
        return surf_vec[surf_index];
      }
    }
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
    double zdriftlength = (cluster->getLocalY() -  _tpc_tzero + _sampa_tzero_bias - 
			   crossing_tzero_correction) * _drift_velocity;  // cm
    double zloc = _max_driftlength/2.0 - zdriftlength;         // local z relative to surface center (for north side):
    unsigned int side = TpcDefs::getSide(key);
    if (side == 0)
    {
      zloc = -zloc;
    }
    local(0) = cluster->getLocalX();
    local(1) = zloc;
  }
  else
  {
    local(0) = cluster->getLocalX();
    local(1) = cluster->getLocalY();
  }

  return local;
}

/*!
 * \file TpcClusterMover.cc
 * \brief moves distortion corrected clusters back to their TPC surface
 * \author Tony Frawley, April 2022
 */

#include "TpcClusterMover.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <iostream>

namespace
{
  [[maybe_unused]] std::ostream& operator<<(std::ostream& out, const Acts::Vector3& v)
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }
}  // namespace

TpcClusterMover::TpcClusterMover()
  : inner_tpc_spacing((mid_tpc_min_radius - inner_tpc_min_radius) / 16.0)
  , mid_tpc_spacing((outer_tpc_min_radius - mid_tpc_min_radius) / 16.0)
  , outer_tpc_spacing((outer_tpc_max_radius - outer_tpc_min_radius) / 16.0)
{
  // initialize layer radii

  for (int i = 0; i < 16; ++i)
  {
    layer_radius[i] = inner_tpc_min_radius + (double) i * inner_tpc_spacing + 0.5 * inner_tpc_spacing;
  }
  for (int i = 0; i < 16; ++i)
  {
    layer_radius[i + 16] = mid_tpc_min_radius + (double) i * mid_tpc_spacing + 0.5 * mid_tpc_spacing;
  }
  for (int i = 0; i < 16; ++i)
  {
    layer_radius[i + 32] = outer_tpc_min_radius + (double) i * outer_tpc_spacing + 0.5 * outer_tpc_spacing;
  }
}

void TpcClusterMover::initialize_geometry(PHG4TpcGeomContainer* cellgeo, ActsGeometry *tGeometry, PHCompositeNode* topNode)
{
  _topNode = topNode;
  
  if (_verbosity > 0)
  {
    std::cout << "TpcClusterMover: Getting ActsGeometry, and getting layer radii for Tpc from cell geometry object" << std::endl;
  }

  if(!tGeometry || !cellgeo)
    {
      std::cout << PHWHERE << " Failed to get ActsGeometry or TPC cell geometry, cannot continue - quit!" << std::endl;
	exit(1);
    }
  
  _tGeometry = tGeometry;
  
  int layer = 0;
  PHG4TpcGeomContainer::ConstRange layerrange = cellgeo->get_begin_end();
  for (PHG4TpcGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    layer_radius[layer] = layeriter->second->get_radius();
    layer++;
  }
}

//____________________________________________________________________________..
std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> TpcClusterMover::processTrack(const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>>& global_in)
{
  // Get the global positions of the TPC clusters for this track, already corrected for distortions, and move them to the surfaces
  // The input object contains all clusters for the track in world coordinates
  // The surface radii are in global coordinates. Needed because we are dealing with aligned surfaces

  auto *cluster_map = findNode::getClass<TrkrClusterContainer>(_topNode, "TRKR_CLUSTER_SEED");
  if(!cluster_map)
    {
      cluster_map = findNode::getClass<TrkrClusterContainer>(_topNode, "TRKR_CLUSTER");
      if(!cluster_map)
	{
	  std::cout << PHWHERE << " Did not find TRKR_CLUSTER_SEED or TRKR_CLUSTER nodes, quit." << std::endl;
	  exit(1);
	}
    }

  auto surfMaps = _tGeometry->maps();
    
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_moved;

  std::vector<Acts::Vector3> tpc_global_vec;
  std::vector<TrkrDefs::cluskey> tpc_cluskey_vec;

  for (const auto& [ckey, global] : global_in)
  {
    const auto trkrid = TrkrDefs::getTrkrId(ckey);
    if (trkrid == TrkrDefs::tpcId)
    {
      tpc_cluskey_vec.push_back(ckey);
      tpc_global_vec.push_back(global);
    }
    else
    {
      // si clusters stay where they are
      global_moved.emplace_back(ckey, global);
    }
  }

  // need at least 3 clusters to fit a circle
  if (tpc_global_vec.size() < 3)
  {
    if (_verbosity > 0)
    {
      std::cout << "  -- skip this tpc track, not enough clusters: " << tpc_global_vec.size() << std::endl;
    }
    return global_in;
  }

  std::vector<float> fitpars = TrackFitUtils::fitClusters(tpc_global_vec,  tpc_cluskey_vec, 0);
  if(fitpars.size() == 0)
    {
      if(_verbosity > 1)
	{
	  std::cout << PHWHERE << "Warning: fit failed, return input positions. " << std::endl;
	}
      return global_in;
    }  
  // Now we need to move each TPC cluster associated with this track to the readout surface radius
  for (unsigned int i = 0; i < tpc_global_vec.size(); ++i)
  {
    TrkrDefs::cluskey cluskey = tpc_cluskey_vec[i];
    Acts::Vector3 global = tpc_global_vec[i];

    // get target surface radius in global coordinates
    auto *cluster = cluster_map->findCluster(cluskey);
    if(!cluster) { continue; }
    auto surface = surfMaps.getSurface(cluskey, cluster);
    if(!surface) { continue; }
    
    // Find the intersection of the helix fitted to the global cluster positions with
    // the surface associated with this cluster key
    TrkrDefs::subsurfkey sskey = cluster->getSubSurfKey();
    Acts::Vector3 global_new= global;
    TrkrDefs::subsurfkey new_subsurfkey = sskey;
    bool ret = get_moved_position(cluskey, cluster, fitpars, global, global_new, new_subsurfkey);
    if(!ret)
    {
      global_moved.emplace_back(cluskey, global); 
      if(_verbosity > 1)
	{
	  std::cout << PHWHERE << "Warning: get_moved_position failed, use input position. " << std::endl;
	}
      continue;
    }
    
    int iter = 0;
    while(new_subsurfkey != sskey)
      {
	iter++;
	if(iter > 2) { break; }

	// surface changed, cluster subsurface has been updated, redo with new surface
	sskey = new_subsurfkey;
	ret = get_moved_position(cluskey, cluster, fitpars, global, global_new, new_subsurfkey);	  
	if(_verbosity > 2)
	  {
	    if(new_subsurfkey != sskey)
	      {
		std::cout << PHWHERE << "Warning: subsurfkey changed on iteration " << iter << " from "
			  << sskey << " to " << new_subsurfkey << std::endl;
	      }
	    else
	      {
		std::cout << PHWHERE << "Good: subsurfkey unchanged on iteration " << iter << std::endl;
	      }
	  }
	
      }

    if(_verbosity > 2)
      {
	std::cout << "    iterations " << iter << " clusterkey " << cluskey << " subsurfkey in " << sskey << " new " << new_subsurfkey << std::endl;
	std::cout << "        global_in " << global[0] << "  " << global[1] << "  " << global[2]  << std::endl;
	std::cout << "        global_moved " << global_new[0] << "  " << global_new[1] << "  " << global_new[2]  << std::endl;
      }

    // add the new position and surface to the return object
    global_moved.emplace_back(cluskey, global_new); 
  }
  
  return global_moved;
}

bool TpcClusterMover::get_moved_position(TrkrDefs::cluskey cluskey, TrkrCluster *cluster, std::vector<float> &fitpars, Acts::Vector3 &global, Acts::Vector3 &global_new, TrkrDefs::subsurfkey &new_subsurfkey) const
{
  auto surfMaps = _tGeometry->maps();
  auto surface = surfMaps.getSurface(cluskey, cluster);
      
  Acts::Vector3 surf_intercept = TrackFitUtils::get_helix_surface_intersection(surface,  fitpars, global, _tGeometry);
  if(fitpars.size() == 0) { return false; }
  
  // get circle position at cluster radius
  double cluster_radius = sqrt(global[0] * global[0] + global[1] * global[1]);
  double R = fitpars[0];
  double X0 = fitpars[1];
  double Y0 = fitpars[2];
  double x_start = 0.0;
  double y_start = 0.0;
  int ret = get_circle_circle_intersection(cluster_radius, R, X0, Y0, global[0], global[1], x_start, y_start);
  if (ret == Fun4AllReturnCodes::ABORTEVENT)
    {
      return false;  // skip to next cluster
    }
  // z projection is unique
  double A = fitpars[3];
  double B = fitpars[4];
  double z_start = B + A * cluster_radius;
  
  Acts::Vector3 start_point(x_start, y_start, z_start); 
  
  // calculate dx, dy, dz along circle trajectory from cluster radius to surface radius
  double xnew = global[0] - (start_point.x() - surf_intercept.x());
  double ynew = global[1] - (start_point.y() - surf_intercept.y());
  double znew = global[2] - (start_point.z() - surf_intercept.z());
  
  global_new(0) = xnew;
  global_new(1) = ynew;
  global_new(2) = znew;
  
  bool update_sskey = true;
  if(update_sskey)
    {
      unsigned int sskey = cluster->getSubSurfKey();
      
      TrkrDefs::hitsetkey hkey= TrkrDefs:: getHitSetKeyFromClusKey(cluskey);
      auto new_surf = _tGeometry->get_tpc_surface_from_coords(hkey, global_new, new_subsurfkey);
      if( new_surf && (new_subsurfkey != sskey) )
	{
	  cluster->setSubSurfKey(new_subsurfkey);
	}
    }
  
  return true;
}

int TpcClusterMover::get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double& x, double& y) const
{
  // finds the intersection of the fitted circle with the cylinder having radius = target_radius
  const auto [xplus, yplus, xminus, yminus] = TrackFitUtils::circle_circle_intersection(target_radius, R, X0, Y0);

  // We only need to check xplus for failure, skip this TPC cluster in that case
  if (std::isnan(xplus))
  {
    {
      if (_verbosity > 1)
      {
        std::cout << " circle/circle intersection calculation failed, skip this cluster" << std::endl;
        std::cout << " target_radius " << target_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
      }
    }
    return Fun4AllReturnCodes::ABORTEVENT;  // skip to next cluster
  }

  // we can figure out which solution is correct based on the cluster position in the TPC
  if (fabs(xclus - xplus) < 5.0 && fabs(yclus - yplus) < 5.0)  // 5 cm, large and arbitrary
  {
    x = xplus;
    y = yplus;
  }
  else
  {
    x = xminus;
    y = yminus;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

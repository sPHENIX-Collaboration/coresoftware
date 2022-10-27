/*!
 * \file TpcClusterMover.cc
 * \brief moves distortion corrected clusters back to their TPC surface
 * \author Tony Frawley, April 2022 
 */

#include "TpcClusterMover.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <trackbase/TrackFitUtils.h>

#include <cmath>
#include <iostream>
#include <climits>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

TpcClusterMover::TpcClusterMover()
{
  // initialize layer radii
  inner_tpc_spacing = (mid_tpc_min_radius - inner_tpc_min_radius) / 16.0;
  mid_tpc_spacing = (outer_tpc_min_radius - mid_tpc_min_radius) / 16.0;
  outer_tpc_spacing = (outer_tpc_max_radius - outer_tpc_min_radius) / 16.0;
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i] = inner_tpc_min_radius + (double) i * inner_tpc_spacing + 0.5 * inner_tpc_spacing;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+16] = mid_tpc_min_radius + (double) i * mid_tpc_spacing + 0.5 * mid_tpc_spacing;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+32] = outer_tpc_min_radius + (double) i * outer_tpc_spacing  +  0.5 * outer_tpc_spacing;
    }
}

void TpcClusterMover::initialize_geometry(PHG4TpcCylinderGeomContainer* cellgeo)
{
  
  int layer=0;
  PHG4TpcCylinderGeomContainer::ConstRange layerrange = cellgeo->get_begin_end();
  for (PHG4TpcCylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    layer_radius[layer] = layeriter->second->get_radius();
    layer++;
  }

}

//____________________________________________________________________________..
std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> TpcClusterMover::processTrack(std::vector<std::pair<TrkrDefs::cluskey,Acts::Vector3>> global_in )
{

  // Get the global positions of the TPC clusters for this track, already corrected for distortions, and move them to the surfaces
  // The input object contains all clusters for the track
    
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_moved;

  std::vector<Acts::Vector3> tpc_global_vec;
  std::vector<TrkrDefs::cluskey> tpc_cluskey_vec;

  for(unsigned int i=0; i< global_in.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = global_in[i].first;
      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
      if(trkrid == TrkrDefs::tpcId)
	{
	  tpc_global_vec.push_back(global_in[i].second);
	  tpc_cluskey_vec.push_back(global_in[i].first);
	}
      else
	{
	  // si clusters stay where they are
	  global_moved.push_back(std::make_pair(cluskey, global_in[i].second));
	}
    }
 
  // need at least 3 clusters to fit a circle
  if(tpc_global_vec.size() < 3)
    {
      if(_verbosity > 0)
	{ std::cout << "  -- skip this tpc track, not enough clusters: " << tpc_global_vec.size() << std::endl; }
      return global_in;
    }
	  
  // fit a circle to the TPC clusters
  const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin( tpc_global_vec );
  
  // get the straight line representing the z trajectory in the form of z vs radius
  const auto [A, B] = TrackFitUtils::line_fit( tpc_global_vec );
  
  // Now we need to move each TPC cluster associated with this track to the readout layer radius
  for(unsigned int i=0; i< tpc_global_vec.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = tpc_cluskey_vec[i];
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      Acts::Vector3 global = tpc_global_vec[i];
      
      // get circle position at target surface radius 
      double target_radius = layer_radius[layer-7];
      int ret = get_circle_circle_intersection(target_radius, R, X0, Y0, global[0], global[1], _x_proj, _y_proj);
      if(ret == Fun4AllReturnCodes::ABORTEVENT) continue;  // skip to next cluster
      // z projection is unique
      _z_proj = B + A * target_radius;
      
      // get circle position at cluster radius	  
      double cluster_radius = sqrt(global[0] * global[0] + global[1] * global[1]);
      ret = get_circle_circle_intersection(cluster_radius, R, X0, Y0, global[0], global[1], _x_start, _y_start);
      if(ret == Fun4AllReturnCodes::ABORTEVENT) continue;  // skip to next cluster
      // z projection is unique
      _z_start = B + A * cluster_radius;
      
      // calculate dx, dy, dz along circle trajectory from cluster radius to surface radius
      double xnew = global[0] - (_x_start - _x_proj);
      double ynew = global[1] - (_y_start - _y_proj);
      double znew = global[2] - (_z_start - _z_proj);
      
      // now move the cluster to the surface radius
      // we keep the cluster key fixed, change the surface if necessary
      
      Acts::Vector3 global_new(xnew, ynew, znew);
      
      // add the new position and surface to the return object
      global_moved.push_back(std::make_pair(cluskey, global_new));
     }

  return global_moved;
}

int TpcClusterMover::get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double &x, double &y)
{
  // finds the intersection of the fitted circle with the cylinder having radius = target_radius
  const auto [xplus, yplus, xminus, yminus] = TrackFitUtils::circle_circle_intersection(target_radius, R, X0, Y0 );
   
  // We only need to check xplus for failure, skip this TPC cluster in that case
   if(std::isnan(xplus)) 
     {
	 {
	   if(_verbosity > 1)
	     {
	       std::cout << " circle/circle intersection calculation failed, skip this cluster" << std::endl;
	       std::cout << " target_radius " << target_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
	     }
	 }
       return Fun4AllReturnCodes::ABORTEVENT;  // skip to next cluster
     }
   
   // we can figure out which solution is correct based on the cluster position in the TPC
   if(fabs(xclus - xplus) < 5.0 && fabs(yclus - yplus) < 5.0)  // 5 cm, large and arbitrary 
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

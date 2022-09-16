#include "PHTpcClusterMover.h"

#include "PHTpcClusterMover.h"   

/// Tracking includes

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrClusterv3.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/ActsTransformations.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

#include <TF1.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
PHTpcClusterMover::PHTpcClusterMover(const std::string &name)
  : SubsysReco(name)
{}

//____________________________________________________________________________..
int PHTpcClusterMover::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  for (int layer=7; layer<7+48; layer++){
    PHG4TpcCylinderGeom* GeoLayer = _tpc_geom_container->GetLayerCellGeom(layer);
    std::cout << "PHTpcClusterMover:: layer = " << layer << " layer_radius " << GeoLayer->get_radius() << std::endl;  
    layer_radius[layer-7] = GeoLayer->get_radius();
  }

  //// initialize layer radii
  //inner_tpc_spacing = (mid_tpc_min_radius - inner_tpc_min_radius) / 16.0;
  //mid_tpc_spacing = (outer_tpc_min_radius - mid_tpc_min_radius) / 16.0;
  //outer_tpc_spacing = (outer_tpc_max_radius - outer_tpc_min_radius) / 16.0;
  //for(int i=0; i < 16; ++i)
  //  {
  //    layer_radius[i] = inner_tpc_min_radius + (double) i * inner_tpc_spacing + 0.5 * inner_tpc_spacing;
  //    if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i] << std::endl;
  //  }
  //for(int i=0; i < 16; ++i)
  //  {
  //    layer_radius[i+16] = mid_tpc_min_radius + (double) i * mid_tpc_spacing + 0.5 * mid_tpc_spacing;
  //    if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i+16] << std::endl;
  //  }
  //for(int i=0; i < 16; ++i)
  //  {
  //    layer_radius[i+32] = outer_tpc_min_radius + (double) i * outer_tpc_spacing  +  0.5 * outer_tpc_spacing;
  //     if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i+32] << std::endl;
  //  }

  return ret;
}

//____________________________________________________________________________..
int PHTpcClusterMover::process_event(PHCompositeNode */*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;

  // loop over the tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      _track = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout << std::endl
	    << __LINE__
	    << ": Processing track itrack: " << phtrk_iter->first
	    << ": nhits: " << _track-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _track->get_phi()
		    << std::endl;
	}

      // Get the TPC clusters for this track and correct them for distortions
      std::vector<Acts::Vector3> globalClusterPositions;
      std::map<TrkrDefs::cluskey, Acts::Vector3> tpc_clusters;

      for (auto key_iter = _track->begin_cluster_keys();
	   key_iter != _track->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  // non Tpc clusters are copied unchanged to the new map if they are present
	  // This is needed for truth seeding case, where the tracks already have silicon clusters 
	  if(trkrId != TrkrDefs::tpcId) 
	    {
        
        // check if clusters has not been inserted already
        if( _corrected_cluster_map->findCluster(cluster_key) ) continue;
        
	      // get cluster from original map
	      auto cluster = _cluster_map->findCluster(cluster_key);	
	      if( !cluster ) continue;
	      
        // create a copy
        auto newclus = new TrkrClusterv3;
        newclus->CopyFrom( cluster );
        
        // insert in corrected map
        _corrected_cluster_map->addClusterSpecifyKey(cluster_key, newclus);
	      continue;      
	    }
	  
	  // get the cluster in 3D coordinates
	  auto tpc_clus =  _cluster_map->findCluster(cluster_key);
	  auto global = _tGeometry->getGlobalPosition(cluster_key, tpc_clus);

	  // check if TPC distortion correction are in place and apply
	  if(Verbosity() > 2)  std::cout << "  layer " << layer << " distorted cluster position: " << global[0] << "  " << global[1] << "  " << global[2];
	  if( _dcc ) global = _distortionCorrection.get_corrected_position( global, _dcc ); 
	  if(Verbosity() > 2) std::cout << "   corrected cluster position: " << global[0] << "  " << global[1] << "  " << global[2] << std::endl;

	  // Store the corrected 3D cluster positions
	  globalClusterPositions.push_back(global);
	  tpc_clusters.insert(std::make_pair(cluster_key, global));

	}

      // need at least 3 clusters to fit a circle
      if(globalClusterPositions.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc track, not enough clusters: " << globalClusterPositions.size() << std::endl; 
	  continue;  // skip to the next TPC track
	}

  // fit a circle to the clusters
  const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin( globalClusterPositions );
  if(Verbosity() > 10) 
  { std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl; }

  // toss tracks for which the fitted circle could not have come from the vertex
  //if(R < 30.0) continue;
  
  // get the straight line representing the z trajectory in the form of z vs radius
  const auto [A, B] = TrackFitUtils::line_fit( globalClusterPositions );
  if(Verbosity() > 10)  
  { std::cout << " Fitted line has A " << A << " B " << B << std::endl; }

  // Now we need to move each cluster associated with this track to the readout layer radius
  for( const auto& [cluskey, global]:tpc_clusters )
  {
    const unsigned int layer = TrkrDefs::getLayer(cluskey);

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
	  // write the new cluster position local coordinates on the surface

	  Acts::Vector3 global_new(xnew, ynew, znew);
	  
	  TrkrDefs::subsurfkey subsurfkey;
	  TrkrDefs::hitsetkey tpcHitSetKey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  Surface surface = _tGeometry->get_tpc_surface_from_coords(
            tpcHitSetKey,
	    global_new,
	    subsurfkey);
	
	  if(!surface)
	    {
	      /// If the surface can't be found, we can't track with it. So 
	      /// just continue and don't modify the cluster to the container
	      std::cout << PHWHERE << "Failed to find surface for cluster " << cluskey << std::endl;
	      continue;
	    }

	  // get the original cluster
	  TrkrCluster *cluster =  _cluster_map->findCluster(cluskey);	

	  // put the corrected cluster in the new cluster map

	  // ghost tracks can have repeat clusters, so check if cluster already moved
	  if(_corrected_cluster_map->findCluster(cluskey)) continue;

    // create new cluster
    auto newclus = new TrkrClusterv3;

    // copy from source
    newclus->CopyFrom( cluster );

    // assign subsurface key
    newclus->setSubSurfKey(subsurfkey);

	  // get local coordinates
    Acts::Vector3 normal = surface->normal(_tGeometry->geometry().getGeoContext());
    auto local = surface->globalToLocal(_tGeometry->geometry().getGeoContext(),
					global * Acts::UnitConstants::cm,
					normal);

	  Acts::Vector2 localPos;
	  if(local.ok())
	    {
	      localPos = local.value() / Acts::UnitConstants::cm;
	    }
	  else
	    {
	      /// otherwise take the manual calculation
	      Acts::Vector3 center = surface->center(_tGeometry->geometry().getGeoContext())/Acts::UnitConstants::cm;
	      double clusRadius = sqrt(xnew * xnew + ynew * ynew);
	      double clusphi = atan2(ynew, xnew);
	      double rClusPhi = clusRadius * clusphi;
	      double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
	      double surfPhiCenter = atan2(center[1], center[0]);
	      double surfRphiCenter = surfPhiCenter * surfRadius;
	      double surfZCenter = center[2];
	      
	      localPos(0) = rClusPhi - surfRphiCenter;
	      localPos(1) = znew - surfZCenter; 
	    }
	  
	  if(Verbosity() > 4)
	    {
	      std::cout << "*** cluster_radius " << cluster_radius << " cluster x,y,z: " << global[0] << "  " << global[1] << "  " << global[2] << std::endl;
	      std::cout << "    projection_radius " << target_radius << " proj x,y,z: " << _x_proj << "  " << _y_proj << "  " << _z_proj << std::endl; 
	      std::cout << "    traj_start_radius " << cluster_radius << " start x,y,z: "<< _x_start << "  " << _y_start << "  " << _z_start << std::endl; 
	      std::cout << "    moved_clus_radius " << target_radius << " final x,y,z: "<< xnew << "  " << ynew << "  " << znew << std::endl; 
	    }

    // assign to new cluster
    newclus->setLocalX(localPos(0));
	  newclus->setLocalY(localPos(1));

    // insert in map
    _corrected_cluster_map->addClusterSpecifyKey(cluskey, newclus);
	}

      // For normal reconstruction, the silicon clusters  for this track will be copied over after the matching is done
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcClusterMover::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTpcClusterMover::GetNodes(PHCompositeNode* topNode)
{
  _tpc_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_tpc_geom_container)
  {
    std::cout << PHWHERE << " ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( _dcc )
    { 
      std::cout << "PHTpcClusterMover:   found TPC distortion correction container" << std::endl; 
    }
      
  // create the node for distortion corrected clusters, if it does not already exist
  _corrected_cluster_map  = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if(!_corrected_cluster_map)
    {
      std::cout << "Creating node CORRECTED_TRKR_CLUSTER" << std::endl;

      PHNodeIterator iter(topNode);

      // Looking for the DST node
      PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
      if (!dstNode)
	{
	  std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}      
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode("TRKR");
	  dstNode->addNode(DetNode);
	}
      
      _corrected_cluster_map = new TrkrClusterContainerv4;
      PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(_corrected_cluster_map, "CORRECTED_TRKR_CLUSTER", "PHObject");
      DetNode->addNode(TrkrClusterContainerNode);
    }    
  else
    {
      _corrected_cluster_map->Reset();
    }          


  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcClusterMover::get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double &x, double &y)
{
  // finds the intersection of the fitted circle with the cylinder having radius = target_radius
  const auto [xplus, yplus, xminus, yminus] = TrackFitUtils::circle_circle_intersection(target_radius, R, X0, Y0 );
   
   // We only need to check xplus for failure, skip this TPC cluster in that case
   if(std::isnan(xplus)) 
     {
       if(Verbosity() > 0)
	 {
	   std::cout << " circle/circle intersection calculation failed, skip this cluster" << std::endl;
	   std::cout << " target_radius " << target_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
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

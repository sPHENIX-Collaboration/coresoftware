#include "PHTpcClusterMover.h"

#include "PHTpcClusterMover.h"   

/// Tracking includes

#include <trackbase/TrkrClusterv3.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/ActsTransformations.h>


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
{

}

//____________________________________________________________________________..
PHTpcClusterMover::~PHTpcClusterMover()
{

}

//____________________________________________________________________________..
int PHTpcClusterMover::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // initialize layer radii
  inner_tpc_spacing = (mid_tpc_min_radius - inner_tpc_min_radius) / 16.0;
  mid_tpc_spacing = (outer_tpc_min_radius - mid_tpc_min_radius) / 16.0;
  outer_tpc_spacing = (outer_tpc_max_radius - outer_tpc_min_radius) / 16.0;
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i] = inner_tpc_min_radius + (double) i * inner_tpc_spacing + 0.5 * inner_tpc_spacing;
      if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i] << std::endl;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+16] = mid_tpc_min_radius + (double) i * mid_tpc_spacing + 0.5 * mid_tpc_spacing;
      if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i+16] << std::endl;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+32] = outer_tpc_min_radius + (double) i * outer_tpc_spacing  +  0.5 * outer_tpc_spacing;
       if(Verbosity() > 4) std::cout << " i " << i << " layer_radius " << layer_radius[i+32] << std::endl;
    }

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
        _corrected_cluster_map->addCluster(newclus);
	      continue;      
	    }
	  
	  // get the cluster in 3D coordinates
	  auto tpc_clus =  _cluster_map->findCluster(cluster_key);
	  auto global = _transformer.getGlobalPosition(tpc_clus,
						      _surfmaps,
						      _tGeometry);

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
      double R = 0;
      double X0 = 0;
      double Y0 = 0;
      CircleFitByTaubin(globalClusterPositions, R, X0, Y0);
      if(Verbosity() > 10) 
	std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      //if(R < 30.0) continue;

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit(globalClusterPositions, A, B);
      if(Verbosity() > 10) 
	std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Now we need to move each cluster associated with this track to the readout layer radius
      for (auto clus_iter = tpc_clusters.begin();
	   clus_iter != tpc_clusters.end(); 
	   ++clus_iter)
	{
	  TrkrDefs::cluskey cluskey = clus_iter->first;
	  unsigned int layer = TrkrDefs::getLayer(cluskey);
	  Acts::Vector3 global = clus_iter->second;

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
	  Surface surface = get_tpc_surface_from_coords(tpcHitSetKey,
							global_new,
							_surfmaps,
							_tGeometry,
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
	  Acts::Vector3 normal = surface->normal(_tGeometry->geoContext);
	  auto local = surface->globalToLocal(_tGeometry->geoContext,
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
	      Acts::Vector3 center = surface->center(_tGeometry->geoContext)/Acts::UnitConstants::cm;
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
    _corrected_cluster_map->addCluster(newclus);
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

  _surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  if(!_surfmaps)
    {
      std::cout << PHWHERE << "Error, can't find acts surface maps" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
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
      
      _corrected_cluster_map = new TrkrClusterContainerv3;
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
  double xplus = 0;
  double yplus = 0; 
   double xminus = 0;
   double yminus = 0;
   
   circle_circle_intersection(target_radius, R, X0, Y0, xplus, yplus, xminus, yminus);
   
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
void PHTpcClusterMover::CircleFitByTaubin (std::vector<Acts::Vector3> clusters, double &R, double &X0, double &Y0)
/*  
      Circle fit to a given set of data points (in 2D)
      This is an algebraic fit, due to Taubin, based on the journal article
      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With 
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
     It works well whether data points are sampled along an entire circle or along a small arc. 
     It still has a small bias and its statistical accuracy is slightly lower than that of the geometric fit (minimizing geometric distances),
     It provides a very good initial guess for a subsequent geometric fit. 
       Nikolai Chernov  (September 2012)
*/
{
  int iter,IterMAX=99;
  
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double x,y;
  double DET,Xcenter,Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;
  for(unsigned int iclus = 0; iclus < clusters.size(); ++iclus)
    {
      if(Verbosity() > 3)  std::cout << "    add cluster with x " << clusters[iclus][0] << " and y " << clusters[iclus][1] << std::endl;
      meanX += clusters[iclus][0];
      meanY += clusters[iclus][1];
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<clusters.size(); i++)
    {
      double Xi = clusters[i][0] - meanX;   //  centered x-coordinates
      double Yi = clusters[i][1] - meanY;   //  centered y-coordinates
      double Zi = Xi*Xi + Yi*Yi;
      
      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;
  
  //  computing coefficients of the characteristic polynomial
  
  Mz = Mxx + Myy;
  Cov_xy = Mxx*Myy - Mxy*Mxy;
  Var_z = Mzz - Mz*Mz;
  A3 = 4*Mz;
  A2 = -3*Mz*Mz - Mzz;
  A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;
  
  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0  
  //    (it is guaranteed to converge to the right root)
  
  for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x*(A22 + A33*x);
      double xnew = x - y/Dy;
      if ((xnew == x)||(!std::isfinite(xnew))) break;
      double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);
}

void PHTpcClusterMover::circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus)
{
  // r1 is radius of sPHENIX layer
  // r2, x2 and y2 are parameters of circle fitted to TPC clusters
  // the solutions are xplus, xminus, yplus, yminus

  // The intersection of two circles occurs when
  // (x-x1)^2 + (y-y1)^2 = r1^2,  / (x-x2)^2 + (y-y2)^2 = r2^2
  // Here we assume that circle 1 is an sPHENIX layer centered on x1=y1=0, and circle 2 is arbitrary
  //  x^2 +y^2 = r1^2,   (x-x2)^2 + (y-y2)^2 = r2^2
  // expand the equations and subtract to eliminate the x^2 and y^2 terms, gives the radical line connecting the intersection points
  // iy = - (2*x2*x - D) / 2*y2, 
  // then substitute for y in equation of circle 1

  double D = r1*r1 - r2*r2 + x2*x2 + y2*y2;
  double a = 1.0 + (x2*x2) / (y2*y2);
  double b = - D * x2/( y2*y2);
  double c = D*D / (4.0*y2*y2) - r1*r1;

  xplus = (-b + sqrt(b*b - 4.0* a * c) ) / (2.0 * a);
  xminus = (-b - sqrt(b*b - 4.0* a * c) ) / (2.0 * a);

  // both values of x are valid
  // but for each of those values, there are two possible y values on circle 1
  // but only one of those falls on the radical line:

  yplus = - (2*x2*xplus - D) / (2.0*y2); 
  yminus = -(2*x2*xminus - D) / (2.0*y2);

}

void  PHTpcClusterMover::line_fit(std::vector<Acts::Vector3> clusters, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius

   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<clusters.size(); ++i)
    {
      double z = clusters[i][2];
      double r = sqrt(pow(clusters[i][0],2) + pow(clusters[i][1], 2));

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(clusters.size()*xysum-xsum*ysum)/(clusters.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*clusters.size()-xsum*xsum);            //calculate intercept

    return;
}   

Surface PHTpcClusterMover::get_tpc_surface_from_coords(TrkrDefs::hitsetkey hitsetkey,
						       Acts::Vector3 world,
						       ActsSurfaceMaps *surfMaps,
						       ActsTrackingGeometry *tGeometry,
						       TrkrDefs::subsurfkey& subsurfkey)
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  std::map<unsigned int, std::vector<Surface>>::iterator mapIter;
  mapIter = surfMaps->tpcSurfaceMap.find(layer);
  
  if(mapIter == surfMaps->tpcSurfaceMap.end())
    {
      std::cout << PHWHERE 
		<< "Error: hitsetkey not found in clusterSurfaceMap, hitsetkey = "
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
    nsurf += surf_vec.size()/2;

  Surface this_surf = surf_vec[nsurf];
      
  auto vec3d = this_surf->center(tGeometry->geoContext);
  std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
  double surf_z = surf_center[2];
  double surf_phi = atan2(surf_center[1], surf_center[0]);
  double surfStepPhi = tGeometry->tpcSurfStepPhi;
  double surfStepZ = tGeometry->tpcSurfStepZ;

  if( (world_phi > surf_phi - surfStepPhi / 2.0 && world_phi < surf_phi + surfStepPhi / 2.0 ) &&
      (world_z > surf_z - surfStepZ / 2.0 && world_z < surf_z + surfStepZ / 2.0) )	
    {
      if(Verbosity() > 2)
	std::cout <<  "     got it:  surf_phi " << surf_phi << " surf_z " << surf_z 
		  << " surfStepPhi/2 " << surfStepPhi/2.0 << " surfStepZ/2 " << surfStepZ/2.0  
		  << " world_phi " << world_phi << " world_z " << world_z 
		  << " rounded_nsurf "<< rounded_nsurf << " surf_index " << nsurf
		  << std::endl;        
      
      surf_index = nsurf;
      subsurfkey = nsurf;
    }    
  else
    {
      std::cout << PHWHERE 
		<< "Error: TPC surface index not defined, skipping cluster!" 
		<< std::endl;
      std::cout << "     coordinates: " << world[0] << "  " << world[1] << "  " << world[2] 
		<< " radius " << sqrt(world[0]*world[0]+world[1]*world[1]) << std::endl;
      std::cout << "     world_phi " << world_phi << " world_z " << world_z << std::endl;
      std::cout << "     surf coords: " << surf_center[0] << "  " << surf_center[1] << "  " << surf_center[2] << std::endl;
      std::cout << "     surf_phi " << surf_phi << " surf_z " << surf_z << std::endl; 
      std::cout << " number of surfaces " << surf_vec.size() << std::endl;
      return nullptr;
    }
  
  return surf_vec[surf_index];
  
}

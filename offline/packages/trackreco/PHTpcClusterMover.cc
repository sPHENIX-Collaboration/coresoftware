#include "PHTpcClusterMover.h"

#include "PHTpcClusterMover.h"   

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

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
  double inner_tpc_spacing = (mid_tpc_min_radius - inner_tpc_min_radius) / 16.0;
  double mid_tpc_spacing = (outer_tpc_min_radius - mid_tpc_min_radius) / 16.0;
  double outer_tpc_spacing = (outer_tpc_max_radius - outer_tpc_min_radius) / 16.0;
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i] = inner_tpc_min_radius + (double) i * inner_tpc_spacing + 0.5 * inner_tpc_spacing;
      if(Verbosity() > 0) std::cout << " i " << i << " layer_radius " << layer_radius[i] << std::endl;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+16] = mid_tpc_min_radius + (double) i * mid_tpc_spacing + 0.5 * mid_tpc_spacing;
      if(Verbosity() > 0) std::cout << " i " << i << " layer_radius " << layer_radius[i+16] << std::endl;
    }
  for(int i=0; i < 16; ++i)
    {
      layer_radius[i+32] = outer_tpc_min_radius + (double) i * outer_tpc_spacing  +  0.5 * outer_tpc_spacing;
       if(Verbosity() > 0) std::cout << " i " << i << " layer_radius " << layer_radius[i+32] << std::endl;
    }

  return ret;
}

//____________________________________________________________________________..
int PHTpcClusterMover::process_event(PHCompositeNode *topNode)
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

      // Get the TPC clusters for this track
      std::map<unsigned int, TrkrCluster*> tpc_clusters;
      std::vector<TrkrCluster*> clusters;

      for (SvtxTrack::ConstClusterKeyIter key_iter = _track->begin_cluster_keys();
	   key_iter != _track->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(trkrId != TrkrDefs::tpcId) continue;

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  tpc_clusters.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 10) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << tpc_clusters.size() << std::endl;
	}

      // need at least 3 clusters to fit a circle
      if(tpc_clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc track, not enough clusters " << std::endl; 
	  continue;  // skip to the next TPC track
	}

      // fit a circle to the clusters
      double R = 0;
      double X0 = 0;
      double Y0 = 0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) 
	std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 30.0) continue;

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit(clusters, A, B);
      if(Verbosity() > 10) 
	std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Now we need to move each cluster associated with this track to the readout layer radius
      for (auto clus_iter = tpc_clusters.begin();
	   clus_iter != tpc_clusters.end(); 
	   ++clus_iter)
	{
	  unsigned int layer = clus_iter->first;
	  TrkrCluster *cluster = clus_iter->second;
	  
	  double target_radius = layer_radius[layer-7];
	  
	  // finds the intersection of the fitted circle with the micromegas layer
	  double xplus = 0;
	  double yplus = 0; 
	  double xminus = 0;
	  double yminus = 0;
	  circle_circle_intersection(target_radius, R, X0, Y0, xplus, yplus, xminus, yminus);
	  
	// We only need to check xplus for failure, skip this TPC track in that case
	if(std::isnan(xplus)) 
	  {
	    if(Verbosity() > 0)
	      {
		std::cout << " circle/circle intersection calculation failed, skip this cluster" << std::endl;
		std::cout << " target_radius " << target_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
	      }
	    continue;  // skip to next cluster
	  }

	// we can figure out which solution is correct based on the cluster position in the TPC
	double xclus = cluster->getX();
	double yclus = cluster->getY();
	if(fabs(xclus - xplus) < 5.0 && fabs(yclus - yplus) < 5.0)
	  {
	    _x_proj = xplus;
	    _y_proj = yplus;
	  }
	else
	  {
	    _x_proj = xminus;
	    _y_proj = yminus;
	  }

	// z projection is unique
	_z_proj = B + A * target_radius;

	if(Verbosity() > 10)
	  {
	    std::cout << "target_radius " << target_radius << " cluster x,y,z:  old " << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ() 
		      << " new " << _x_proj << "  " << _y_proj << "  " << _z_proj << std::endl; 
	  }

	// now move the cluster to the projected position
	cluster->setX(_x_proj);
	cluster->setY(_y_proj);
	cluster->setZ(_z_proj);
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcClusterMover::End(PHCompositeNode *topNode)
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
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTpcClusterMover::CircleFitByTaubin (std::vector<TrkrCluster*> clusters, double &R, double &X0, double &Y0)
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
      meanX += clusters[iclus]->getX();
      meanY += clusters[iclus]->getY();
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<clusters.size(); i++)
    {
      double Xi = clusters[i]->getX() - meanX;   //  centered x-coordinates
      double Yi = clusters[i]->getY() - meanY;   //  centered y-coordinates
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
      if ((xnew == x)||(!isfinite(xnew))) break;
      double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (abs(ynew)>=abs(y))  break;
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

void  PHTpcClusterMover::line_fit(std::vector<TrkrCluster*> clusters, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<clusters.size(); ++i)
    {
      double z = clusters[i]->getZ();
      double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(clusters.size()*xysum-xsum*ysum)/(clusters.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*clusters.size()-xsum*xsum);            //calculate intercept

   if(Verbosity() > 10)
     {
       for (unsigned int i=0;i<clusters.size(); ++i)
	 {
	   double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
	   double z_fit = a * r + b;                    //to calculate y(fitted) at given x points
	   std::cout << " r " << r << " z " << clusters[i]->getZ() << " z_fit " << z_fit << std::endl; 
	 } 
     }

    return;
}   

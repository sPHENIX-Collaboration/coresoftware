/*!
 * \file TpcClusterMover.cc
 * \brief moves distortion corrected clusters back to their TPC surface
 * \author Tony Frawley, April 2022 
 */

#include "TpcClusterMover.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <cmath>
#include <iostream>
#include <climits>

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

//____________________________________________________________________________..
std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> TpcClusterMover::processTrack(std::vector<std::pair<TrkrDefs::cluskey,Acts::Vector3>> global_in )
{

  // Get the TPC clusters for this track and correct them for distortions
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
      std::cout << "  -- skip this tpc track, not enough clusters: " << tpc_global_vec.size() << std::endl; 
      return global_in;
    }
	  
  // fit a circle to the TPC clusters
  double R = 0;
  double X0 = 0;
  double Y0 = 0;
  CircleFitByTaubin(tpc_global_vec, R, X0, Y0);

  //std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;
  
  // get the straight line representing the z trajectory in the form of z vs radius
  double A = 0; double B = 0;
  line_fit(tpc_global_vec, A, B);

  //std::cout << " Fitted line has A " << A << " B " << B << std::endl;
  
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

void TpcClusterMover::CircleFitByTaubin (std::vector<Acts::Vector3> clusters, double &R, double &X0, double &Y0)
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

void  TpcClusterMover::line_fit(std::vector<Acts::Vector3> clusters, double &a, double &b)
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

void TpcClusterMover::circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus)
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

int TpcClusterMover::get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double &x, double &y)
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

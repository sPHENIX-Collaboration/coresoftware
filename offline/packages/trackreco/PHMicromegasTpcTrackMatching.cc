#include "PHMicromegasTpcTrackMatching.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>

#include "AssocInfoContainer.h"

#include <TF1.h>

using namespace std;

//____________________________________________________________________________..
PHMicromegasTpcTrackMatching::PHMicromegasTpcTrackMatching(const std::string &name):
 PHTrackPropagating(name)
{

}

//____________________________________________________________________________..
PHMicromegasTpcTrackMatching::~PHMicromegasTpcTrackMatching()
{

}

//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::Setup(PHCompositeNode *topNode)
{

  std::cout << std::endl << PHWHERE 
	    << " rphi_search_win_1 " << _rphi_search_win_1
	    << " z_search_win_1 " << _z_search_win_1
	    << " rphi_search_win_2 " << _rphi_search_win_2
	    << " z_search_win_2 " << _z_search_win_2
	    << endl;

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::Process()
{
  // _track_map contains the TPC seed track stubs
  // We will add the micromegas cluster to the TPC tracks already on the node tree
  // We will have to expand the number of tracks whenever we find multiple matches to the silicon

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size() << endl;

 // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;

  // loop over the original TPC tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;
      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout << std::endl
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << endl;
	}

      // Get the outermost TPC clusters for this tracklet
      std::map<unsigned int, TrkrCluster*> outer_clusters;
      std::vector<TrkrCluster*> clusters;

      for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
	   key_iter != _tracklet_tpc->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(layer < _min_tpc_layer) continue;
	  if(layer >= _min_mm_layer) continue;

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  outer_clusters.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 3) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << outer_clusters.size() << std::endl;
	}
      
      if(outer_clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough outer clusters " << std::endl; 
	  continue;
	}

      // fit a circle to the clusters
      double R, X0, Y0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 3) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;


      // make the function to extrapolate to the micromegas layers
      // find the most extreme layers
      TrkrCluster *clus_layers[2];      
      unsigned int max_layer = 0;
      unsigned int min_layer = 100;
      TrkrCluster *clus_min_layer = 0;      
      TrkrCluster *clus_max_layer = 0;

      for(auto clus_iter=outer_clusters.begin(); clus_iter != outer_clusters.end(); ++clus_iter)
	{
	  unsigned int layer = clus_iter->first;	  
	  if(layer > max_layer)
	    { 
	      max_layer = layer;
	      clus_max_layer = clus_iter->second;	      
	    }
	  if(layer < min_layer) 
	    {
	      min_layer = layer;
	      clus_min_layer = clus_iter->second;
	    }
	}

      clus_layers[0] = clus_min_layer;
      clus_layers[1] = clus_max_layer;

      // get the straight line representing the z trajectory in the form of z vs radius
      double zclus1= clus_min_layer->getZ();
      double rclus1 = sqrt(clus_min_layer->getX()*clus_min_layer->getX() + clus_min_layer->getY()*clus_min_layer->getY() );
      double zclus2 = clus_max_layer->getZ();
      double rclus2 = sqrt(clus_max_layer->getX()*clus_max_layer->getX() + clus_max_layer->getY()*clus_max_layer->getY() );

      // define the straight line joining the two clusters
      // make sure the cluster in the inner layer is in [0]
      float xl[2] = {clus_layers[0]->getX(), clus_layers[1]->getX()};
      float yl[2] = {clus_layers[0]->getY(), clus_layers[1]->getY()};
      float zl[2] = {clus_layers[0]->getZ(), clus_layers[1]->getZ()};

      // loop over the micromegas clusters and find any within the search windows
      std::vector<TrkrDefs::cluskey> mm_matches_1;
      std::vector<TrkrDefs::cluskey> mm_matches_2;
      //TrkrClusterContainer::ConstRange mm_clusrange = _cluster_map->getClusters(TrkrDefs::micromegasId);
      TrkrClusterContainer::ConstRange mm_clusrange = _cluster_map->getClusters();
      for(TrkrClusterContainer::ConstIterator clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
	{
	  TrkrDefs::cluskey mm_cluskey = clusiter->first;
	  unsigned int layer = TrkrDefs::getLayer(mm_cluskey);
	  TrkrCluster *mm_clus = clusiter->second;

	  if(layer < _min_mm_layer) continue;

	  if(Verbosity() > 3)
	    {
	      std::cout << " Found Micromegas cluster in layer " << layer  << " radius, x, y, z " 
			<< sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2)) << "  "
			<< mm_clus->getX() << "  " << mm_clus->getY() << "  " << mm_clus->getZ() << std::endl;
	    }

	  double mm_clus_z = mm_clus->getZ();
	  double mm_radius = sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2) );
	  double mm_clus_rphi = mm_radius * atan2(mm_clus->getY(), mm_clus->getX());

	  float x_proj, y_proj, z_proj, rphi_proj,radius_proj;

	  // method to find where fitted circle intersects this layer
	  double xplus = 0;
	  double xminus = 0;
	  double yplus = 0;
	  double yminus = 0;
	  // finds the intersection of the fitted circle with the micromegas layer
	  circle_circle_intersection(mm_radius, R, X0, Y0, xplus, yplus, xminus, yminus);
	  // choose the answer closest to the micromegas cluster in phi
	  double rphi_plus = mm_radius* atan2(yplus, xplus);
	  double rphi_minus = mm_radius * atan2(yminus, xminus);
	  double x_proj_circle;
	  double y_proj_circle;
	  double rphi_proj_circle;
	  if( fabs(rphi_plus - mm_clus_rphi) < fabs(rphi_minus - mm_clus_rphi) )
	    {
	      rphi_proj_circle = rphi_plus; 
	      x_proj_circle = xplus;
	      y_proj_circle = yplus;
	    }
	  else
	    {
	      rphi_proj_circle = rphi_minus; 
	      x_proj_circle = xminus;
	      y_proj_circle = yminus;
	    }

	  // understand why this happens occasionally
	  if(isnan(rphi_proj_circle)) 
	    {
	      std::cout << " rphi_proj_circle is nan. xplus " << xplus << " yplus " << yplus << " xminus " << xminus << " yminus " << yminus << std::endl; 
	      continue;
	    }

	  double radius_proj_circle;
	  radius_proj_circle = sqrt(x_proj_circle*x_proj_circle + y_proj_circle*y_proj_circle);

	  // find where z projects to on the micromegas layer
	  double z_proj_circle = zclus2 + (mm_radius - rclus2) * (zclus2 - zclus1) / (rclus2 - rclus1) ;

	  if(Verbosity() > 3)
	    std::cout << "     circle:  radius_proj " << radius_proj_circle << " x_proj " << x_proj_circle 
		      << " y_proj " << y_proj_circle << " z_proj " << z_proj_circle  
		      << " rphi_proj " << rphi_proj_circle << std::endl;

	  float t = NAN;	  
	  t = line_circle_intersection(xl, yl, zl, mm_radius);
	  if (t > 0)
	    {
	      x_proj = xl[0] + t * (xl[1] - xl[0]);
	      y_proj = yl[0] + t * (yl[1] - yl[0]);
	      z_proj = zl[0] + t * (zl[1] - zl[0]);
	      rphi_proj = atan2(y_proj, x_proj) * mm_radius;	  
	      radius_proj = sqrt(x_proj*x_proj + y_proj * y_proj);

	      if(Verbosity() > 3)
		{
		  std::cout << "     linear:  radius_proj " << radius_proj << " x_proj " << x_proj 
			    << " y_proj " << y_proj << " z_proj " << z_proj << " rphi_proj " << rphi_proj 
			    << std::endl;
		}
	    }
	  else
	    {
	      std::cout << "    line-circle intersetion fit failed, skip this micromegas cluster " << endl;
	    continue;
	    }

	  bool use_circle = true;
	  if(use_circle)
	    {
	      rphi_proj = rphi_proj_circle;
	      z_proj = z_proj_circle;
	    }

	  if(layer == 55)
	    {
	      if(Verbosity() > 3)
		{
		  std::cout << "     test for match in layer " << layer << " _rphi_search_win_1 " << _rphi_search_win_1
			    << " drphi " << rphi_proj - mm_clus_rphi << " dz " << z_proj - mm_clus_z
			    << " _z_search_win_1 " << _z_search_win_1 
			    << std::endl;
		}
	      
	      std::cout << "     deltas " << layer  << " drphi " << rphi_proj - mm_clus_rphi << " dz " << z_proj - mm_clus_z << std::endl;
	      
	      if(fabs(rphi_proj - mm_clus_rphi) < _rphi_search_win_1 && fabs(z_proj - mm_clus_z) < _z_search_win_1)
		{
		  mm_matches_1.push_back(mm_cluskey);
		}
	    }
	  
	  if(layer == 56)
	    {
	      if(Verbosity() > 3)
		{
		  std::cout << "     test for match in layer " << layer << " _rphi_search_win_2 " << _rphi_search_win_2
			    << " drphi " << rphi_proj - mm_clus_rphi << " dz " << z_proj - mm_clus_z
			    << " _z_search_win_2 " << _z_search_win_2 
			    << std::endl;
		}
	      
	      std::cout << "    deltas " << layer  << " drphi " << rphi_proj - mm_clus_rphi << " dz " << z_proj - mm_clus_z << std::endl;
	      
	      if(fabs(z_proj - mm_clus_z) < _z_search_win_2 && fabs(rphi_proj - mm_clus_rphi) < _rphi_search_win_2)
		{
		  mm_matches_2.push_back(mm_cluskey);
		}
	    }	  
	}
      
      // skip no match or multiple match cases for now
      if(mm_matches_1.size() != 1) continue;
      if(mm_matches_2.size() != 1) continue;
      
      if(Verbosity() > 3)
	{
	  cout << "Original TPC tracklet:" << endl;
	  _tracklet_tpc->identify();
	}
      
      // Add the micromegas clusters to the track
      
      if(Verbosity() > 3)
	std::cout << "  insert cluster key " << mm_matches_1[0] << " into tracklet " << _tracklet_tpc->get_id() << std::endl;
      /*
      // don't do this until surfaces are in Acts geometry
      _tracklet_tpc->insert_cluster_key(mm_matches_1[0]);
      _assoc_container->SetClusterTrackAssoc(mm_matches_1[0], _tracklet_tpc->get_id());
      */
      
      if(Verbosity() > 3)
	std::cout << "  insert cluster key " << mm_matches_2[0] << " into tracklet " << _tracklet_tpc->get_id() << std::endl;
      
      /*
      // don't do this until surfaces are in Acts geometry
      _tracklet_tpc->insert_cluster_key(mm_matches_2[0]);
      _assoc_container->SetClusterTrackAssoc(mm_matches_2[0], _tracklet_tpc->get_id());
      */
      if(Verbosity() > 3)
	_tracklet_tpc->identify();
      
    }
  
  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() << endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHMicromegasTpcTrackMatching::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHMicromegasTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  
  return Fun4AllReturnCodes::EVENT_OK;
}

float  PHMicromegasTpcTrackMatching::line_circle_intersection(float x[], float y[], float z[], float radius)
{
  // parameterize the line in terms of t (distance along the line segment, from 0-1) as
  // x = x0 + t * (x1-x0); y=y0 + t * (y1-y0); z = z0 + t * (z1-z0)
  // parameterize the cylinder (centered at x,y = 0,0) as  x^2 + y^2 = radius^2,   then
  // (x0 + t*(x1-z0))^2 + (y0+t*(y1-y0))^2 = radius^2
  // (x0^2 + y0^2 - radius^2) + (2x0*(x1-x0) + 2y0*(y1-y0))*t +  ((x1-x0)^2 + (y1-y0)^2)*t^2 = 0 = C + B*t + A*t^2
  // quadratic with:  A = (x1-x0)^2+(y1-y0)^2 ;  B = 2x0*(x1-x0) + 2y0*(y1-y0);  C = x0^2 + y0^2 - radius^2
  // solution: t = (-B +/- sqrt(B^2 - 4*A*C)) / (2*A)

  float A = (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]);
  float B = 2.0 * x[0] * (x[1] - x[0]) + 2.0 * y[0] * (y[1] - y[0]);
  float C = x[0] * x[0] + y[0] * y[0] - radius * radius;
  float tup = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
  float tdn = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

  // The limits are 0 and 1, but we allow a little for floating point precision
  float t;
  if (tdn >= 0)
    t = tdn;
  else if (tup > 0)
    t = tup;
  else
  {
    cout << PHWHERE << "   **** Oops! No valid solution for tup or tdn, tdn = " << tdn << " tup = " << tup << endl;
    cout << "   radius " << radius << " rbegin " << sqrt(x[0] * x[0] + y[0] * y[0]) << " rend " << sqrt(x[1] * x[1] + y[1] * y[1]) << endl;
    cout << "   x0 " << x[0] << " x1 " << x[1] << endl;
    cout << "   y0 " << y[0] << " y1 " << y[1] << endl;
    cout << "   z0 " << z[0] << " z1 " << z[1] << endl;
    cout << "   A " << A << " B " << B << " C " << C << endl;

    t = -1;
  }

  return t;
}

void PHMicromegasTpcTrackMatching::CircleFitByTaubin (std::vector<TrkrCluster*> clusters, double &R, double &X0, double &Y0)
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
  int i,iter,IterMAX=99;
  
  double Xi,Yi,Zi;
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double Dy,xnew,x,ynew,y;
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
  
  for (i=0; i<clusters.size(); i++)
    {
      Xi = clusters[i]->getX() - meanX;   //  centered x-coordinates
      Yi = clusters[i]->getY() - meanY;   //  centered y-coordinates
      Zi = Xi*Xi + Yi*Yi;
      
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
      Dy = A1 + x*(A22 + A33*x);
      xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
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

void PHMicromegasTpcTrackMatching::circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus)
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


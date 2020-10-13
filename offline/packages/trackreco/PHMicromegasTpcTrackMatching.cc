#include "PHMicromegasTpcTrackMatching.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <micromegas/MicromegasDefs.h>

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
	    << " rphi_search_win inner layer " << _rphi_search_win[0]
	    << " z_search_win inner layer " << _z_search_win[0]
	    << " rphi_search_win outer layer " << _rphi_search_win[1]
	    << " z_search_win outer layer " << _z_search_win[1]
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

  _event++;

  if(Verbosity() > 0)
    cout << PHWHERE << " Event " << _event << " TPC track map size " << _track_map->size() << endl;

  // get the micromegas clusters for this event
  TrkrDefs::TrkrId mm_trkrid = TrkrDefs::micromegasId;; 
  TrkrClusterContainer::ConstRange mm_clusrange = _cluster_map->getClusters(mm_trkrid);

  // for diagnostics only
  std::map<TrkrDefs::cluskey, int> cluster_matches;
  for(TrkrClusterContainer::ConstIterator clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
    {
      TrkrDefs::cluskey mm_cluskey = clusiter->first;
      cluster_matches.insert(make_pair(mm_cluskey, 0));
    }
  
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

	  if(Verbosity() > 10) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << outer_clusters.size() << std::endl;
	}


      // need at least 3 clusters to fit a circle
      if(outer_clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough outer clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      // fit a circle to the clusters
      double R, X0, Y0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 40.0) continue;

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit(clusters, A, B);
      if(Verbosity() > 10) std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Project this TPC tracklet  to the two micromegas layers and store the projections
      bool skip_tracklet = false;
      for(unsigned int imm = 0; imm < _n_mm_layers; ++imm)
      {
	// method to find where fitted circle intersects this layer

	double xplus = 0;
	double xminus = 0;
	double yplus = 0;
	double yminus = 0;
	// finds the intersection of the fitted circle with the micromegas layer
	circle_circle_intersection(	_mm_layer_radius[imm], R, X0, Y0, xplus, yplus, xminus, yminus);
	
	// We only need to check xplus for failure, skip this TPC track in that case
	if(std::isnan(xplus)) 
	  {
	    if(Verbosity() > 10)
	      {
		std::cout << " circle/circle intersection calculation failed, skip this case" << std::endl;
		std::cout << " mm_radius " << _mm_layer_radius[imm] << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
	      }
	    skip_tracklet = true;
	  }
	if(skip_tracklet == true)
	  continue;   // skips to the next layer

	// we can figure out which solution is correct based on the last cluster position in the TPC
	unsigned int nlast = clusters.size() -1;
	double last_clus_phi = atan2(clusters[nlast]->getY(), clusters[nlast]->getX());
	double plus_phi = atan2(yplus, xplus);
	double minus_phi = atan2(yminus, xminus);
	if(fabs(last_clus_phi - plus_phi) < fabs(last_clus_phi - minus_phi))
	  {
	    _rphi_proj[imm] = plus_phi * _mm_layer_radius[imm]; 
	    _x_proj[imm] = xplus;
	    _y_proj[imm] = yplus;
	  }
	else
	  {
	    _rphi_proj[imm] = minus_phi * _mm_layer_radius[imm]; 
	    _x_proj[imm] = xminus;
	    _y_proj[imm] = yminus;
	  }

	// z projection is unique
	_z_proj[imm] = B + A * _mm_layer_radius[imm] ;
      }
      
      if(skip_tracklet == true)
	continue;   // skips to the next TPC tracklet
      
      // loop over the micromegas clusters and find any within the search windows
      std::vector<TrkrDefs::cluskey> mm_matches[2];
       for(TrkrClusterContainer::ConstIterator clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
	{
	  TrkrDefs::cluskey mm_cluskey = clusiter->first;
	  unsigned int layer = TrkrDefs::getLayer(mm_cluskey);
	  TrkrCluster *mm_clus = clusiter->second;

	  unsigned int imm;
	  if(layer == _min_mm_layer) 
	    imm = 0;
	  else
	    imm = 1;
	  
	  if(Verbosity() > 3)
	    {
	      std::cout << " Found Micromegas cluster in layer " << layer  << " cluskey " << mm_cluskey << " radius, x, y, z, phi, rphi " 
			<< sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2)) << "  "
			<< mm_clus->getX() << "  " << mm_clus->getY() << "  " << mm_clus->getZ() << "  " 
			<< atan2(mm_clus->getY(), mm_clus->getX()) << "  "
			<<  sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2))  * atan2(mm_clus->getY(), mm_clus->getX()) << std::endl;
	    }

	  double mm_clus_z = mm_clus->getZ();
	  double mm_radius = sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2) );
	  double mm_clus_rphi = mm_radius * atan2(mm_clus->getY(), mm_clus->getX());

	  double radius_proj = sqrt(_x_proj[imm]*_x_proj[imm] + _y_proj[imm]*_y_proj[imm]);
	  
	  if(Verbosity() > 3)
	    {
	      std::cout << "   tracklet " << _tracklet_tpc->get_id() << " test for match in layer " << layer << " _rphi_search_win_1 " << _rphi_search_win[imm]
			<< " phi_proj " << _rphi_proj[imm] / radius_proj << " drphi " << _rphi_proj[imm] - mm_clus_rphi << " _z_proj " << _z_proj[imm] 
			<< " dz " << _z_proj[imm] - mm_clus_z
			<< " _z_search_win " << _z_search_win[imm] 
			<< std::endl;
	    }
	  
	  if(fabs(_rphi_proj[imm] - mm_clus_rphi) < _rphi_search_win[imm] && fabs(_z_proj[imm] - mm_clus_z) < _z_search_win[imm])
	    {
	      mm_matches[imm].push_back(mm_cluskey);
	      cluster_matches.find(mm_cluskey)->second++;

	      if(Verbosity() > 3)
		std::cout << "     radius_proj " << radius_proj << " _x_proj " << _x_proj 
			  << " _y_proj " << _y_proj << " _z_proj " << _z_proj[imm]  
			  << " _rphi_proj " << _rphi_proj[imm] << std::endl;
	      

	      // prints out a line that can be grep-ed from the output file to feed to a display macro
	      if( _test_search_windows )
		std::cout << "     deltas " << layer  << " drphi " << _rphi_proj[imm] - mm_clus_rphi << " dz " << _z_proj[imm] - mm_clus_z 
			  << " mm_clus_rphi " << mm_clus_rphi << " mm_clus_z " << mm_clus_z << " match " << mm_matches[imm].size()  << std::endl;
	    }
	}

       // We need to modify the Micromegas cluster position for the unmeasured coordinates so they have the projected position, instead of the middle of the tile
       // this is so that the cluster will be associated in PHActsSourceLinks with the surface that Acts will project the track to
       // This kludge is only necessary until we implement the Micromegas tiles properly in Geant
       
       for(unsigned int imm = 0; imm < _n_mm_layers; ++imm)
	 {      
	   for(unsigned int imatch = 0; imatch < mm_matches[imm].size(); ++imatch)
	     {
	       TrkrCluster *cluster = _cluster_map->findCluster(mm_matches[imm][imatch]);
	       
	       // update the coordinate that is not measured
	       const auto segmentationType(MicromegasDefs::getSegmentationType(mm_matches[imm][imatch]));
	       switch( segmentationType )
		 {
		 case MicromegasDefs::SegmentationType::SEGMENTATION_PHI: 	      
		   cluster->setZ(_z_proj[imm]);
		   break;
		 case MicromegasDefs::SegmentationType::SEGMENTATION_Z: 
		   const auto radius = std::sqrt( pow(cluster->getX(), 2) + pow(cluster->getY(), 2) );
		   cluster->setX(radius*std::cos(_rphi_proj[imm] / radius));
		   cluster->setY(radius*std::sin(_rphi_proj[imm] / radius));
		   break;
		 }
	       
	       if(Verbosity() > 3)
		 std::cout << "       imm " << imm << " imatch " <<  imatch << " updated Micromegas cluster  " << mm_matches[imm][imatch] << " has radius, x, y, z " 
			   << sqrt(pow(cluster->getX(), 2) + pow(cluster->getY(), 2)) << "  "
			   << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ() << std::endl;
	     }
	 }
       
       // keep multiple-matches to TPC tracklets and let the tracxker sort it out
       // but if there are no matches we are done with this track
       if(mm_matches[0].size() == 0 && mm_matches[1].size() == 0) continue;

       if(Verbosity() > 3)
	 {
	   cout << "Original TPC tracklet:" << endl;
	   _tracklet_tpc->identify();
	 }
       
       // Add the micromegas clusters to the track
       for(unsigned int imm = 0; imm < _n_mm_layers; ++imm)
	 {      	  
	   for(unsigned int imatch = 0; imatch < mm_matches[imm].size(); ++imatch)
	     {
	       if(Verbosity() > 3) 
		 std::cout << "   inserting Micromegas cluster with key " << mm_matches[imm][imatch] << std::endl;

	       _tracklet_tpc->insert_cluster_key(mm_matches[imm][imatch]);
	       _assoc_container->SetClusterTrackAssoc(mm_matches[imm][imatch], _tracklet_tpc->get_id());
	     }
	 }
      
       if(Verbosity() > 3)
	 _tracklet_tpc->identify();
    }

  /*
  // temporary diagnostics
  for(TrkrClusterContainer::ConstIterator clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
    {
      TrkrDefs::cluskey mm_cluskey = clusiter->first;
      unsigned int layer = TrkrDefs::getLayer(mm_cluskey);
      TrkrCluster *mm_clus = clusiter->second;
      
      int got_it = cluster_matches.find(mm_cluskey)->second;
      TrkrDefs::hitsetkey this_hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(mm_cluskey);
      unsigned int this_tileid = MicromegasDefs::getTileId(this_hitsetkey);

      if(got_it == 0)      
	std::cout << " got_it " << got_it << " micromegas cluster in layer " << layer  << " cluskey " << mm_cluskey << " tileid " << this_tileid << " radius, x, y, z " 
		  << sqrt(pow(mm_clus->getX(), 2) + pow(mm_clus->getY(), 2)) << "  "
		  << mm_clus->getX() << "  " << mm_clus->getY() << "  " << mm_clus->getZ() << std::endl;
    }
  */
      
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

void  PHMicromegasTpcTrackMatching::line_fit(std::vector<TrkrCluster*> clusters, double &a, double &b)
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

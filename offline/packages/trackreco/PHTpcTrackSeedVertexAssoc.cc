#include "PHTpcTrackSeedVertexAssoc.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>
using namespace std;

//____________________________________________________________________________..
PHTpcTrackSeedVertexAssoc::PHTpcTrackSeedVertexAssoc(const std::string &name):
 PHTrackPropagating(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHTpcTrackSeedVertexAssoc::PHTpcTrackSeedVertexAssoc(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHTpcTrackSeedVertexAssoc::~PHTpcTrackSeedVertexAssoc()
{

}

//____________________________________________________________________________..
int PHTpcTrackSeedVertexAssoc::Setup(PHCompositeNode *topNode)
{
  // put these in the output file
  cout << PHWHERE << " p0 " << _par0 << " p1 " << _par1 << " p2 " 
       << _par2 << " Search windows: phi " << _phi_search_win << " eta " 
       << _eta_search_win << endl;

  // corrects the PHTpcTracker phi bias
  fdphi = new TF1("f1", "[0] + [1]/x^[2]");
  fdphi->SetParameter(0, _par0);
  fdphi->SetParameter(1, _par1);
  fdphi->SetParameter(2, _par2);

  // corrects the space charge distortion phi bias
  if(!_is_ca_seeder)
    {
      // PHTpcTracker correction is opposite in sign
      // and different in magnitude - why?
      _parsc0 *= -1.0 * 0.7;
      _parsc1 *= -1.0 * 0.7;
    }
  fscdphi = new TF1("f2","[0] + [1]*x^2");
  fscdphi->SetParameter(0, _parsc0 * _collision_rate / _reference_collision_rate);
  fscdphi->SetParameter(1, _parsc1 * _collision_rate / _reference_collision_rate);

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHTpcTrackSeedVertexAssoc::Process()
{
  // _track_map contains the TPC seed track stubs
  // We want to associate these TPC track seeds with a collision vertex
  // Then we add the collision vertex position as the track seed position

  // All we need is to project the TPC clusters in Z to the beam line.


  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size()  << endl;
  /*
 // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;
  */

  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      /*
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;
      */      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << ": eta: " << _tracklet_tpc->get_eta()
	    << endl;
	}

      // get the tpc track seed cluster positions in z and r

      // Get the outermost TPC clusters for this tracklet
      std::map<unsigned int, TrkrCluster*> tpc_clusters_map;
      std::vector<TrkrCluster*> clusters;
      
      for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
	   key_iter != _tracklet_tpc->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(layer < _min_tpc_layer) continue;
	  if(layer > _max_tpc_layer) continue;

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  tpc_clusters_map.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 5) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " clusters.size() " << tpc_clusters_map.size() << std::endl;
	}


      // need at least 3 clusters to fit a circle
      if(tpc_clusters_map.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      /*
      // fit a circle to the clusters
      double R, X0, Y0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 40.0) continue;
      */

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit_clusters(clusters, A, B);
      if(Verbosity() > 5) std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Project this TPC tracklet  to the beam line and store the projections
      //bool skip_tracklet = false;
      // z projection is unique
      _z_proj = B;
      
      // Find the nearest collision vertex
      // if find multiple close matches, duplicate the track?

      int trackVertexId = 9999;
      double dz = 9999.;	  
      for(SvtxVertexMap::Iter viter = _vertex_map->begin();
	  viter != _vertex_map->end();
	  ++viter)
	{
	  auto vertexKey = viter->first;
	  auto vertex = viter->second;
	  if(Verbosity() > 100)
	    vertex->identify();

	  const double vertexZ = vertex->get_z();
	  
	  if(fabs(_z_proj - vertexZ) < dz )
	    {
	      dz = fabs(_z_proj - vertexZ);
	      trackVertexId = vertexKey;
	    }	  
	}  // end loop over collision vertices

      _tracklet_tpc->set_vertex_id(trackVertexId);
      auto vertex = _vertex_map->find(trackVertexId)->second;
      vertex->insert_track(phtrk_iter->first);

      // set the track position to the vertex position
      _tracklet_tpc->set_x(vertex->get_x());
      _tracklet_tpc->set_y(vertex->get_y());
      _tracklet_tpc->set_z(vertex->get_z());

      if(Verbosity() > 1)
	{
	  std::cout << "TPC seed track " << phtrk_iter->first << " matched to vertex " << trackVertexId << endl; 
	} 

      // Finished association of track with vertex
      // Now we modify the track parameters

      // Repeat the z line fit including the vertex position, get theta, update pz
      std::vector<std::pair<double, double>> points;
      double r_vertex = sqrt(vertex->get_x()*vertex->get_x() + vertex->get_y()*vertex->get_y());
      double z_vertex = vertex->get_z();
      points.push_back(make_pair(r_vertex, z_vertex));
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double z = clusters[i]->getZ();
	  double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
	  
	  points.push_back(make_pair(r,z));
	}
      
      line_fit(points, A, B);
      if(Verbosity() > 5) 
	std::cout << " Fitted line including vertex has A " << A << " B " << B << std::endl;      

      // extract the track theta
      double track_angle = atan(A);  // referenced to 90 degrees

      // make circle fit including vertex as point
      std::vector<std::pair<double, double>> cpoints;
      double x_vertex = vertex->get_x();
      double y_vertex = vertex->get_y();
      cpoints.push_back(std::make_pair(x_vertex, y_vertex));
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double x = clusters[i]->getX();
	  double y = clusters[i]->getY();	  
	  cpoints.push_back(make_pair(x, y));
	}
      double R, X0, Y0;
      CircleFitByTaubin(cpoints, R, X0, Y0);
      if(Verbosity() > 5) 
	std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      double pt_track = _tracklet_tpc->get_pt();
      bool use_circle_pt = false;    // true for testing only, normally use track seed pT
      if(use_circle_pt)
	{
	  //  get the pT from radius of circle as a check 
	  double Bfield = 1.4;  // T
	  double pt_new = 0.3 * Bfield * R * 0.01;  // convert cm to m 
	  std::cout << " pT from circle of radius R = " << R << " in field of " << Bfield << " Tesla is " << pt_new << " seed pT is " << _tracklet_tpc->get_pt()  << std::endl; 
	  //pt_track = pt_new;
	}

      // We want the angle of the tangent relative to the positive x axis
      // start with the angle of the radial line from vertex to circle center
      double dx = X0 - x_vertex;
      double dy = Y0 - y_vertex;
      double phi= atan2(dy,dx);
      //std::cout << "x_vertex " << x_vertex << " y_vertex " << y_vertex << " X0 " << X0 << " Y0 " << Y0 << " angle " << phi * 180 / 3.14159 << std::endl; 
      // convert to the angle of the tangent to the circle
      // we need to know if the track proceeds clockwise or CCW around the circle
      double dx0 = cpoints[0].first - X0;
      double dy0 = cpoints[0].second - Y0;
      double phi0 = atan2(dy0, dx0);
      double dx1 = cpoints[1].first - X0;
      double dy1 = cpoints[1].second - Y0;
      double phi1 = atan2(dy1, dx1);
      double dphi = phi1 - phi0;

      // need to deal with the switch from -pi to +pi at phi = 180 degrees
      // final phi - initial phi must be < 180 degrees for it to be a valid track
      if(dphi > M_PI) dphi -= 2.0 * M_PI;
      if(dphi < - M_PI) dphi += M_PI;

      if(Verbosity() > 5) 
	{
	  int charge = _tracklet_tpc->get_charge();   // needed for diagnostic output only
	  std::cout << " charge " << charge << " phi0 " << phi0*180.0 / M_PI << " phi1 " << phi1*180.0 / M_PI << " dphi " << dphi*180.0 / M_PI << std::endl;
	}

      // whether we add or subtract 90 degrees depends on the track propagation direction determined above
      if(dphi < 0)
	phi += M_PI / 2.0;  
      else
	phi -= M_PI / 2.0;  
      if(Verbosity() > 5) 
	std::cout << " input track phi " << _tracklet_tpc->get_phi() * 180.0 / M_PI << " new phi " << phi * 180 / M_PI << std::endl;  

      // get the updated values of px, py, pz from the pT and the angles found here
      double px_new = pt_track * cos(phi);
      double py_new = pt_track * sin(phi);
      double ptrack_new = pt_track / cos(track_angle);
      double pz_new = ptrack_new * sin(track_angle);

      if(Verbosity() > 5)
	std::cout << " input track mom " << _tracklet_tpc->get_p() << " new mom " << ptrack_new
		  << " px in " << _tracklet_tpc->get_px()  << " px " << px_new 
		  << " py in " << _tracklet_tpc->get_py() << " py " << py_new 
		  << " pz in " << _tracklet_tpc->get_pz() << " pz " << pz_new 
		  << " eta in " <<  _tracklet_tpc->get_eta() << " phi in " <<  _tracklet_tpc->get_phi() * 180.0 / M_PI
		  << " track angle " << track_angle * 180.0 / M_PI 
		  << std::endl;

      // update track on node tree
      _tracklet_tpc->set_px(px_new);
      _tracklet_tpc->set_py(py_new);
      _tracklet_tpc->set_pz(pz_new);

      if(Verbosity() > 5)
	std::cout << " new mom " <<  _tracklet_tpc->get_p() <<  "  new eta " <<  _tracklet_tpc->get_eta() << " new phi " << _tracklet_tpc->get_phi() * 180.0 / M_PI << std::endl;
      
    }  // end loop over TPC track seeds
  
  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() << endl;
  
  if (Verbosity() >= 1)
    cout << "PHTpcTrackSeedVertexAssoc::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcTrackSeedVertexAssoc::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTpcTrackSeedVertexAssoc::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------
  /*
  _track_map_silicon = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxSiliconTrackMap");
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxSiliconTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */  
  return Fun4AllReturnCodes::EVENT_OK;
}


void  PHTpcTrackSeedVertexAssoc::line_fit(std::vector<std::pair<double,double>> points, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<points.size(); ++i)
    {
      //double z = clusters[i]->getZ();
      //double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
      double r = points[i].first;
      double z = points[i].second;

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(points.size()*xysum-xsum*ysum)/(points.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*points.size()-xsum*xsum);            //calculate intercept

   if(Verbosity() > 10)
     {
       for (unsigned int i=0;i<points.size(); ++i)
	 {
	   double r = points[i].first;
	   double z_fit = a * r + b;                    //to calculate z(fitted) at given r points
	   std::cout << " r " << r << " z " << points[i].second << " z_fit " << z_fit << std::endl; 
	 } 
     }

    return;
}   

void  PHTpcTrackSeedVertexAssoc::line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b)
{
  std::vector<std::pair<double,double>> points;
  
   for (unsigned int i=0; i<clusters.size(); ++i)
     {
       double z = clusters[i]->getZ();
       double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));

       points.push_back(make_pair(r,z));
     }

   line_fit(points, a, b);

    return;
}

void PHTpcTrackSeedVertexAssoc::CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0)
/*  
      Circle fit to a given set of data points (in 2D)
      This is an algebraic fit, due to Taubin, based on the journal article
      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With 
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
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
  for(unsigned int i = 0; i < points.size(); ++i)
    {
      meanX += points[i].first;
      meanY += points[i].second;
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<points.size(); i++)
    {
      double Xi = points[i].first - meanX;   //  centered x-coordinates
      double Yi = points[i].second - meanY;   //  centered y-coordinates
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

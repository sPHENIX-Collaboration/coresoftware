#include "PHTpcTrackSeedCircleFit.h"

#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v3.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/ActsTransformations.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

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
#include <cmath>                              // for fabs, std::sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

namespace
{
   
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // radius
  template<class T> T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  void line_fit(std::vector<std::pair<double,double>> points, double &a, double &b)
  {
    // copied from: https://www.bragitoff.com
    // we want to fit z vs radius
    
    //variables for sums/sigma of xi,yi,xi^2,xiyi etc
    double xsum=0,x2sum=0,ysum=0,xysum=0;        
    for( const auto& [r,z]:points )
    {      
      xsum+=r;                        //calculate sigma(xi)
      ysum+=z;                        //calculate sigma(yi)
      x2sum+=square(r);                //calculate sigma(x^2i)
      xysum+=r*z;                    //calculate sigma(xi*yi)
    }
    a=(points.size()*xysum-xsum*ysum)/(points.size()*x2sum-square(xsum)); //calculate slope
    b=(x2sum*ysum-xsum*xysum)/(points.size()*x2sum-square(xsum));         //calculate intercept
  
    return;
  }
  
  void line_fit_clusters(std::vector<Acts::Vector3>& globPos, double &a, double &b)
  {
    std::vector<std::pair<double,double>> points;
    std::transform( globPos.begin(), globPos.end(), std::back_inserter( points ), []( const Acts::Vector3& pos ) 
      { return std::make_pair( get_r( pos(0), pos(1) ), pos(2)); } );

    line_fit(points, a, b);
  }
  
  void CircleFitByTaubin (std::vector<Acts::Vector3> points, double &R, double &X0, double &Y0)
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
    for(const auto& position : points)
    {
      meanX += position(0);
      meanY += position(1);
      weight++;
    }
    meanX /= weight;
    meanY /= weight;
    
    //     computing moments 
    
    Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
    
    for (const auto& position : points)
    {
      double Xi = position(0) - meanX;   //  centered x-coordinates
      double Yi = position(1) - meanY;   //  centered y-coordinates
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
    
    for (x=0.,y=A0,iter=0; iter<IterMAX; ++iter)  // usually, 4-6 iterations are enough
    {
      const double Dy = A1 + x*(A22 + A33*x);
      const double xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      const double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (std::abs(ynew)>=std::abs(y))  break;
      x = xnew;  y = ynew;
    }
    
    //  computing parameters of the fitting circle
    
    DET = x*x - x*Mz + Cov_xy;
    Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
    Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
    
    //  assembling the output
    
    X0 = Xcenter + meanX;
    Y0 = Ycenter + meanY;
    R = std::sqrt(square(Xcenter) + square(Ycenter) + Mz);
  }
  
//   std::vector<double> GetCircleClusterResiduals(std::vector<std::pair<double,double>> points, double R, double X0, double Y0)
//   {
//     std::vector<double> residues;
//     std::transform( points.begin(), points.end(), std::back_inserter( residues ), 
//       [R,X0,Y0]( const std::pair<double,double>& point )
//       {
//         // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
//         const auto& x = point.first;
//         const auto& y = point.second;
//         return std::sqrt( square(x-X0) + square(y-Y0))-R;  
//       } );
//     return residues;  
//   }
  
//   std::vector<double> GetLineClusterResiduals(std::vector<std::pair<double,double>> points, double A, double B)
//   {
//     std::vector<double> residues;
//     std::transform( points.begin(), points.end(), std::back_inserter( residues ), 
//       [A,B]( const std::pair<double,double>& point )
//       {
//         // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
//         const auto& r = point.first;
//         const auto& z = point.second;
//         const double a = -A;
//         const double b = 1.0;
//         const double c = -B;
//         return std::abs(a*r+b*z+c)/std::sqrt(square(a)+square(b));      
//       });
//     return residues;  
//   }
  
  void findRoot(
    const double R, const double X0,
    const double Y0, double& x,
    double& y)
  {
    /**
    * We need to determine the closest point on the circle to the origin
    * since we can't assume that the track originates from the origin
    * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
    * minimize d = std::sqrt((0-x)^2+(0-y)^2), the distance between the 
    * origin and some (currently, unknown) point on the circle x,y.
    * 
    * Solving the circle eqn for x and substituting into d gives an eqn for
    * y. Taking the derivative and setting equal to 0 gives the following 
    * two solutions. We take the smaller solution as the correct one, as 
    * usually one solution is wildly incorrect (e.g. 1000 cm)
    */
    
    const double miny = (std::sqrt(square(X0)*square(R)*square(Y0) + square(R) 
      * pow(Y0,4)) + square(X0)*Y0 + pow(Y0, 3)) 
      / (square(X0)+square(Y0));
    
    const double miny2 = (-std::sqrt(square(X0)*square(R)*square(Y0) + square(R) 
      * pow(Y0,4)) + square(X0) * Y0 + pow(Y0, 3)) 
      / (square(X0)+square(Y0));

    const double minx = std::sqrt(square(R) - square(miny-Y0)) + X0;
    const double minx2 = -std::sqrt(square(R) - square(miny2-Y0)) + X0;

    /// Figure out which of the two roots is actually closer to the origin
    if(std::abs(minx) < std::abs(minx2))
      x = minx;
    else
      x = minx2;
    
    if(std::abs(miny) < std::abs(miny2))
      y = miny;
    else
      y = miny2;
  }
  
}

//____________________________________________________________________________..
PHTpcTrackSeedCircleFit::PHTpcTrackSeedCircleFit(const std::string &name):
 SubsysReco(name)
{}

//____________________________________________________________________________..
int PHTpcTrackSeedCircleFit::InitRun(PHCompositeNode *topnode)
{ 
  // get relevant nodes
  int ret = GetNodes( topnode );
  if( ret != Fun4AllReturnCodes::EVENT_OK ) return ret;

  return Fun4AllReturnCodes::EVENT_OK; }

//____________________________________________________________________________..
int PHTpcTrackSeedCircleFit::process_event(PHCompositeNode*)
{
  
  // _track_map contains the TPC seed track stubs
  // We want to associate these TPC track seeds with a collision vertex
  // All we need is to project the TPC clusters in Z to the beam line.
  // The TPC track seeds are given a position that is the PCA of the line and circle fit to the beam line

  if(Verbosity() > 0)
    std::cout << PHWHERE << " TPC track map size " << _track_map->size()  << std::endl;

  
  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      const auto& track_key = phtrk_iter->first;
      SvtxTrack* tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << tracklet_tpc-> size_cluster_keys()
	    << ": pT: " << tracklet_tpc->get_pt()
	    << ": phi: " << tracklet_tpc->get_phi()
	    << ": eta: " << tracklet_tpc->get_eta()
	    << std::endl;
	}

      // get the tpc track seed cluster positions in z and r

      // Get the TPC clusters for this tracklet
      std::vector<TrkrCluster*> clusters = getTrackClusters(tracklet_tpc);
      if(Verbosity() > 3) std::cout << " TPC tracklet " << tracklet_tpc->get_id() << " clusters.size " << clusters.size() << std::endl;

      // count TPC layers for this track
      bool reject_track = false;
      std::set<unsigned int> layers;
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  if(!clusters[i])
	    {
	      if(Verbosity() > 0) std::cout << " trackid " << phtrk_iter->first << " no cluster found, skip track" << std::endl;
	      reject_track = true;
	      break;
	    }	      
	  unsigned int layer = TrkrDefs::getLayer(clusters[i]->getClusKey());
	  layers.insert(layer);
	}

      if(reject_track) continue;
    
      unsigned int nlayers = layers.size();
      if(Verbosity() > 2) std::cout << "    TPC layers this track: " << nlayers << std::endl;


      std::vector<Acts::Vector3> globalClusterPositions;
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  const Acts::Vector3 global = getGlobalPosition(clusters.at(i));
	  globalClusterPositions.push_back(global);

	  if(Verbosity() > 3)
	    {
	      ActsTransformations transformer;
	      auto global_before = transformer.getGlobalPosition(clusters.at(i),
								 _surfmaps,
								 _tGeometry);
	      TrkrDefs::cluskey key = clusters.at(i)->getClusKey();
	      std::cout << "CircleFit: Cluster: " << key << " _corrected_clusters " << _are_clusters_corrected << std::endl;
	      std::cout << " Global before: " << global_before[0] << "  " << global_before[1] << "  " << global_before[2] << std::endl;
	      std::cout << " Global after   : " << global[0] << "  " << global[1] << "  " << global[2] << std::endl;
	    }
	}
      
      if(clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough TPC clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; 
      double B = 0;
      line_fit_clusters(globalClusterPositions, A, B);
      if(Verbosity() > 2)
      { std::cout << "PHTpcTrackSeedCircleFit::process_event - track: " << track_key << " A=" << A << " B=" << B << std::endl; }

      // Project this TPC tracklet  to the beam line and store the projections
      const double z_proj = B;
      
      // set the track Z position to the Z dca
      tracklet_tpc->set_z(z_proj);

      // Finished association of track with vertex, now we modify the track parameters

      // extract the track theta
      double track_angle = atan(A);  // referenced to 90 degrees
  
      
      // make circle fit
      double R, X0, Y0;
      CircleFitByTaubin(globalClusterPositions, R, X0, Y0);
      if(Verbosity() > 2) 
      { std::cout << "PHTpcTrackSeedCircleFit::process_event - track: " << track_key << " R=" << R << " X0=" << X0 << " Y0=" << Y0 << std::endl; }

      // set the track x and y positions to the circle PCA
      double dcax, dcay;
      findRoot(R, X0, Y0, dcax, dcay);
      if(std::isnan(dcax) or std::isnan(dcay)) 
	{ continue; }

      tracklet_tpc->set_x(dcax);
      tracklet_tpc->set_y(dcay);

      double pt_track = tracklet_tpc->get_pt();
      if(Verbosity() > 5)
      {
        //  get the pT from radius of circle as a check
        static constexpr double Bfield = 1.4;  // T
        double pt_circle = 0.3 * Bfield * R * 0.01;  // convert cm to m
        std::cout << "PHTpcTrackSeedCircleFit::process_event-"
          << " R: " << R
          << " Bfield: " << Bfield
          << " pt_circle: " << pt_circle << " pt_track: " << pt_track
          << std::endl;
      }
      
      // We want the angle of the tangent relative to the positive x axis
      // start with the angle of the radial line from vertex to circle center
      
      double dx = X0 - dcax;
      double dy = Y0 - dcay;
      double phi= atan2(-dx,dy);
     
      // convert to the angle of the tangent to the circle
      // we need to know if the track proceeds clockwise or CCW around the circle
      double dx0 = globalClusterPositions.at(0)(0) - X0;
      double dy0 = globalClusterPositions.at(0)(1) - Y0;
      double phi0 = atan2(dy0, dx0);
      double dx1 = globalClusterPositions.at(1)(0) - X0;
      double dy1 = globalClusterPositions.at(1)(1) - Y0;
      double phi1 = atan2(dy1, dx1);
      double dphi = phi1 - phi0;

      // need to deal with the switch from -pi to +pi at phi = 180 degrees
      // final phi - initial phi must be < 180 degrees for it to be a valid track
      if(dphi > M_PI) dphi -= 2.0 * M_PI;
      if(dphi < - M_PI) dphi += M_PI;

      if(Verbosity() > 5) 
	{
	  std::cout << " charge " <<  tracklet_tpc->get_charge() << " phi0 " << phi0*180.0 / M_PI << " phi1 " << phi1*180.0 / M_PI 
		    << " dphi " << dphi*180.0 / M_PI << std::endl;
	}

      // whether we add 180 degrees depends on the angle of the bend
      if(dphi < 0)
	{ 
	  phi += M_PI; 
	  if(phi > M_PI)
	    { phi -= 2. * M_PI; }
	}
   
      if(Verbosity() > 5) 
	std::cout << " input track phi " << tracklet_tpc->get_phi()  << " new phi " << phi  << std::endl;  
     
      // get the updated values of px, py, pz from the pT and the angles found here
      double px_new = pt_track * cos(phi);
      double py_new = pt_track * sin(phi);
      double ptrack_new = pt_track / cos(track_angle);
      double pz_new = ptrack_new * sin(track_angle);

      if(Verbosity() > 4)
	std::cout << " input track id " << tracklet_tpc->get_id() 
		  << " input track mom " << tracklet_tpc->get_p() << " new mom " << ptrack_new
		  << " px in " << tracklet_tpc->get_px()  << " px " << px_new 
		  << " py in " << tracklet_tpc->get_py() << " py " << py_new 
		  << " pz in " << tracklet_tpc->get_pz() << " pz " << pz_new 
		  << " eta in " <<  tracklet_tpc->get_eta() << " phi in " <<  tracklet_tpc->get_phi() * 180.0 / M_PI
		  << " track angle " << track_angle * 180.0 / M_PI 
		  << std::endl;

      // update track on node tree
      tracklet_tpc->set_px(px_new);
      tracklet_tpc->set_py(py_new);
      tracklet_tpc->set_pz(pz_new);

      if(Verbosity() > 5)
      {
        std::cout << " new mom " <<  tracklet_tpc->get_p() <<  "  new eta " <<  tracklet_tpc->get_eta()
          << " new phi " << tracklet_tpc->get_phi() * 180.0 / M_PI << std::endl;
      }
    }  // end loop over TPC track seeds

    if(Verbosity() > 0)
      std::cout << " Final track map size " << _track_map->size() << std::endl;

    if (Verbosity() > 0)
      std::cout << "PHTpcTrackSeedCircleFit::process_event(PHCompositeNode *topNode) Leaving process_event" << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcTrackSeedCircleFit::End(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }

int  PHTpcTrackSeedCircleFit::GetNodes(PHCompositeNode* topNode)
{
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

  if(_use_truth_clusters)
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  else
    {
      _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
      if(_cluster_map)
	{
	  std::cout << " using CORRECTED_TRKR_CLUSTER node  "<< std::endl;
	}
      else
	{
	  std::cout << " CORRECTED_TRKR_CLUSTER node not found, using TRKR_CLUSTER " << std::endl;
	  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
	  _are_clusters_corrected = false;
	}
    }
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( _dcc )
  { std::cout << "PHTpcTrackSeedCircleFit::get_Nodes  - found TPC distortion correction container" << std::endl; }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<TrkrCluster*> PHTpcTrackSeedCircleFit::getTrackClusters(SvtxTrack *tracklet_tpc)
{
  std::vector<TrkrCluster*> clusters;
  
  for(SvtxTrack::ConstClusterKeyIter key_iter = tracklet_tpc->begin_cluster_keys(); 
    key_iter != tracklet_tpc->end_cluster_keys();
    ++key_iter)
  {
    const TrkrDefs::cluskey& cluster_key = *key_iter;
    
    // only consider TPC clusters
    if( TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId ) continue;
    
    // get the cluster
    TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);
    
    // insert in list
    clusters.push_back(tpc_clus);
    
    if(Verbosity() > 5) 
    {
      unsigned int layer = TrkrDefs::getLayer(cluster_key);      
      std::cout << "  TPC cluster in layer " << layer << " with local position " << tpc_clus->getLocalX() 
        << "  " << tpc_clus->getLocalY() << " clusters.size() " << clusters.size() << std::endl;
      }
    }
    return clusters;
  }
  
Acts::Vector3 PHTpcTrackSeedCircleFit::getGlobalPosition( TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  ActsTransformations transformer;
  auto globalpos = transformer.getGlobalPosition(cluster,
    _surfmaps,
    _tGeometry);

  // check if TPC distortion correction are in place and apply if clusters are not from the corrected node
  if( !_are_clusters_corrected)
    if(_dcc) { globalpos = _distortionCorrection.get_corrected_position( globalpos, _dcc ); }

  return globalpos;
}


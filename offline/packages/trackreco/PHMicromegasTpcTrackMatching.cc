#include "PHMicromegasTpcTrackMatching.h"

#include "AssocInfoContainer.h"
//#include "PHTrackPropagating.h"     // for PHTrackPropagating

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

/// Tracking includes
#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TVector3.h>

#include <array>
#include <cmath>                              // for sqrt, std::abs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

namespace
{
  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! get radius from x and y
  template<class T>
    inline constexpr T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

}

//____________________________________________________________________________..
PHMicromegasTpcTrackMatching::PHMicromegasTpcTrackMatching(const std::string &name):
 PHTrackPropagating(name)
{}

//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::Setup(PHCompositeNode *topNode)
{

  std::cout << std::endl << PHWHERE 
	    << " rphi_search_win inner layer " << _rphi_search_win[0]
	    << " z_search_win inner layer " << _z_search_win[0]
	    << " rphi_search_win outer layer " << _rphi_search_win[1]
	    << " z_search_win outer layer " << _z_search_win[1]
	    << std::endl;

  // load micromegas geometry
  _geomContainerMicromegas = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );
  if(!_geomContainerMicromegas)
  {
    std::cout << PHWHERE << "Could not find CYLINDERGEOM_MICROMEGAS_FULL." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // ensures there are at least two micromegas layers
  if( _geomContainerMicromegas->get_NLayers() != _n_mm_layers )
  {
    std::cout << PHWHERE << "Inconsistent number of micromegas layers." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // get first micromegas layer
  _min_mm_layer = static_cast<CylinderGeomMicromegas*>(_geomContainerMicromegas->GetFirstLayerGeom())->get_layer();
  
  // fill radii
  for( int i = 0; i < _n_mm_layers; ++i )
  { _mm_layer_radius[i] =  static_cast<CylinderGeomMicromegas*>(_geomContainerMicromegas->GetLayerGeom( _min_mm_layer+i ) )->get_radius(); }
  
  fdrphi = new TF1("fdrphi", "[0] + [1]*std::abs(x)");
  fdrphi->SetParameter(0, _par0 *_collision_rate / _reference_collision_rate);
  fdrphi->SetParameter(1, _par1);

  // evaluation
  _test_windows = true;
  if( _test_windows )
  {
    for( int i = 0; i < _n_mm_layers; ++i )
    { 
      _rphi_residuals[i] = new TH1F( Form( "rphi_%i", i ), Form( "r#Delta#phi layer %i", i ), 100, -1, 1 );
      _rphi_residuals[i]->GetXaxis()->SetTitle( "r#Delta#phi (cm)" );

      _z_residuals[i] = new TH1F( Form( "z_%i", i ), Form( "#Deltaz layer %i", i ), 100, -1, 1 );
      _z_residuals[i]->GetXaxis()->SetTitle( "r#Delta#phi (cm)" );
    }
  }
    
  
  // base class setup
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
    std::cout << PHWHERE << " Event " << _event << " TPC track map size " << _track_map->size() << std::endl;

  // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;

  // loop over the original TPC tracks
  for (auto phtrk_iter = _track_map->begin(); phtrk_iter != _track_map->end(); ++phtrk_iter)
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
        << std::endl;
    }
    
    // Get the outermost TPC clusters for this tracklet
    std::map<unsigned int, TrkrCluster*> outer_clusters;
    std::vector<TrkrCluster*> clusters;
    
    for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys(); key_iter != _tracklet_tpc->end_cluster_keys(); ++key_iter)
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
      {
        std::cout 
          << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
          << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << outer_clusters.size() << std::endl;
      }
    }
    
    // need at least 3 clusters to fit a circle
    if(outer_clusters.size() < 3)
    {
      if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough outer clusters " << std::endl; 
      continue;  // skip to the next TPC tracklet
    }
    
    // fit a circle to the clusters
    double R = 0;
    double X0 = 0;
    double Y0 = 0;
    CircleFitByTaubin(clusters, R, X0, Y0);
    if(Verbosity() > 10) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;
    
    // toss tracks for which the fitted circle could not have come from the vertex
    if(R < 40.0) continue;
    
    // get the straight line representing the z trajectory in the form of z vs radius
    double A = 0; double B = 0;
    line_fit(clusters, A, B);
    if(Verbosity() > 10) std::cout << " Fitted line has A " << A << " B " << B << std::endl;
    
    // store intersection to each micromegas layer
    std::array<TVector3,_n_mm_layers> world_intersection;
    
    // Project this TPC tracklet to the two micromegas layers and store the projections
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
        break;
      }
      
      // we can figure out which solution is correct based on the last cluster position in the TPC
      const double last_clus_phi = atan2(clusters.back()->getY(), clusters.back()->getX());
      const double phi_plus = atan2(yplus, xplus);
      const double phi_minus = atan2(yminus, xminus);
      
      // select the angle that is the closest to last cluster
      const double r = _mm_layer_radius[imm];
      const double z = B + A * _mm_layer_radius[imm];
      
      // store phi, apply coarse space charge corrections in calibration mode
      double phi = std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus) ? phi_plus:phi_minus;
      if(_sc_calib_mode) phi -= fdrphi->Eval(z)/r;
      
      // create intersection point in world coordinates
      world_intersection[imm] = TVector3( r*std::cos(phi), r*std::sin(phi), z );
      
    }   // end loop over Micromegas layers
    
    // skips to the next TPC tracklet
    if(skip_tracklet) continue;  
    
    // loop over the micromegas clusters and find any within the search windows
    std::vector<TrkrDefs::cluskey> mm_matches[2];
    
    const auto mm_hitsetrange = _mm_hitsets->getHitSets(TrkrDefs::micromegasId);
    for( auto hitsetitr = mm_hitsetrange.first; hitsetitr != mm_hitsetrange.second; ++hitsetitr)
    {
      const auto mm_clusrange = _cluster_map->getClusters(hitsetitr->first);
      for(TrkrClusterContainer::ConstIterator clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
      {
        
        const TrkrDefs::cluskey mm_cluskey = clusiter->first;
        const auto layer = TrkrDefs::getLayer(mm_cluskey);
        TrkrCluster *mm_clus = clusiter->second;
        
        const unsigned int imm = (layer == _min_mm_layer) ? 0:1;
        
        const double mm_clus_z = mm_clus->getZ();
        const double mm_radius = get_r( mm_clus->getX(), mm_clus->getY() );
        const double mm_phi = std::atan2( mm_clus->getY(), mm_clus->getX() );
        const double mm_clus_rphi = mm_radius * mm_phi;
        
        if(Verbosity() > 3)
        {
          std::cout << " Found Micromegas cluster in layer " << layer  << " cluskey " << mm_cluskey << " radius, x, y, z, phi, rphi " 
            << mm_radius << "  "
            << mm_clus->getX() << "  " << mm_clus->getY() << "  " << mm_clus->getZ() << "  " 
            << mm_phi << "  "
            << mm_clus_rphi 
            << std::endl;
        }
        
        const double radius_proj = get_r( world_intersection[imm].x(), world_intersection[imm].y() );
        const double phi_proj = std::atan2( world_intersection[imm].y(), world_intersection[imm].x() );
        const double rphi_proj = radius_proj*phi_proj;
        const double z_proj = world_intersection[imm].z();
        
        if(Verbosity() > 3)
        {
          std::cout << "   tracklet " << _tracklet_tpc->get_id() << " test for match in layer " << layer << " _rphi_search_win_1 " << _rphi_search_win[imm]
            << " phi_proj " << phi_proj << " drphi " << rphi_proj - mm_clus_rphi << " _z_proj " << z_proj
            << " dz " << z_proj - mm_clus_z
            << " _z_search_win " << _z_search_win[imm] 
            << std::endl;
        }

        // fill histogram
        if( _test_windows )
        {
          _rphi_residuals[imm]->Fill( rphi_proj - mm_clus_rphi );
          _z_residuals[imm]->Fill( z_proj - mm_clus_z );
        }
        
        // compare to search windows
        if(std::abs( rphi_proj - mm_clus_rphi) < _rphi_search_win[imm] && std::abs(z_proj - mm_clus_z) < _z_search_win[imm] )
        {
          
          // store cluster
          mm_matches[imm].push_back(mm_cluskey);
          
          // prints out a line that can be grep-ed from the output file to feed to a display macro
          if( _test_windows )
          {
            std::cout 
              << "  Try_mms: " << layer  << " drphi " << rphi_proj - mm_clus_rphi  << " dz " << z_proj - mm_clus_z 
              << " mm_clus_rphi " << mm_clus_rphi << " mm_clus_z " << mm_clus_z << " rphi_proj " <<  rphi_proj << " z_proj " << z_proj << std::endl;
          }
          
        }
        
      }
      
    }
    
    // keep multiple-matches to TPC tracklets and let the tracxker sort it out
    // but if there are no matches we are done with this track
    if(mm_matches[0].empty() && mm_matches[1].empty()) continue;
    
    if(Verbosity() > 3)
    {
      std::cout << "Original TPC tracklet:" << std::endl;
      _tracklet_tpc->identify();
    }
    
    // Add the micromegas clusters to the track
    for(unsigned int imm = 0; imm < _n_mm_layers; ++imm)
    {      	  
      for( const auto& key:mm_matches[imm] )
      {
        if(Verbosity() > 3) 
        { std::cout << "   inserting Micromegas cluster with key " << key << std::endl; }
        
        _tracklet_tpc->insert_cluster_key(key);
        _assoc_container->SetClusterTrackAssoc(key, _tracklet_tpc->get_id());
      }
    }
    
    if(Verbosity() > 3)
    { _tracklet_tpc->identify(); }
    
  }
  
  if(Verbosity() > 0)  
  { std::cout << " Final track map size " << _track_map->size() << std::endl; }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________________________
int PHMicromegasTpcTrackMatching::End()
{ 

  if( _test_windows )
  {
    std::unique_ptr<TFile> tfile( TFile::Open( _evaluation_rootfile.c_str(), "RECREATE" ) );
    for( int i = 0; i < _n_mm_layers; ++i )
    { 
      _rphi_residuals[i]->Write();
      _z_residuals[i]->Write();
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK; 
}

//_________________________________________________________________________________________________
int  PHMicromegasTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  
  // micromegas geometry
  _geomContainerMicromegas = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );
  if(!_geomContainerMicromegas)
  {
    std::cout << PHWHERE << "Could not find CYLINDERGEOM_MICROMEGAS_FULL." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // hitsets
  _mm_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!_mm_hitsets)
  {
    std::cout << PHWHERE << "No hitset container on node tree. Bailing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
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
    const double z = clusters[i]->getZ();
    const double r = get_r( clusters[i]->getX(), clusters[i]->getY() );
    
    xsum=xsum+r;                        //calculate sigma(xi)
    ysum=ysum+z;                        //calculate sigma(yi)
    x2sum=x2sum+ square(r);             //calculate sigma(x^2i)
    xysum=xysum+r*z;                    //calculate sigma(xi*yi)
  }
  
  a=(clusters.size()*xysum-xsum*ysum)/(clusters.size()*x2sum-xsum*xsum);            //calculate slope
  b=(x2sum*ysum-xsum*xysum)/(x2sum*clusters.size()-xsum*xsum);            //calculate intercept

  if(Verbosity() > 10)
  {
    for (unsigned int i=0;i<clusters.size(); ++i)
    {
      const double r = get_r( clusters[i]->getX(), clusters[i]->getY() );
      const double z_fit = a * r + b;                    //to calculate y(fitted) at given x points
      std::cout << " r " << r << " z " << clusters[i]->getZ() << " z_fit " << z_fit << std::endl; 
    } 
  }

  return;

}   

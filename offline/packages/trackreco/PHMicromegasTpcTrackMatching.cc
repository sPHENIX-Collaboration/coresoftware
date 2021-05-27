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

#include <TF1.h>
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

  /**
   * Circle fit to a given set of data points (in 2D)
   * This is an algebraic fit, due to Taubin, based on the journal article
   * G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   * Space Curves Defined By Implicit Equations, With
   * Applications To Edge And Range Image Segmentation",
   * IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   * It works well whether data points are sampled along an entire circle or along a small arc.
   * It still has a small bias and its statistical accuracy is slightly lower than that of the geometric fit (minimizing geometric distances),
   * It provides a very good initial guess for a subsequent geometric fit.
   * Nikolai Chernov  (September 2012)
   */
  void CircleFitByTaubin (std::vector<TrkrCluster*> clusters, double &R, double &X0, double &Y0)
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

  // 2D linear fit
  void line_fit(std::vector<TrkrCluster*> clusters, double &a, double &b)
  {
    // copied from: https://www.bragitoff.com
    // we want to fit z vs radius

    //variables for sums/sigma of xi,yi,xi^2,xiyi etc
    double xsum=0,x2sum=0,ysum=0,xysum=0;
    for (unsigned int i=0; i<clusters.size(); ++i)
    {
      const double z = clusters[i]->getZ();
      const double r = get_r( clusters[i]->getX(), clusters[i]->getY() );

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+ square(r);             //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }

    a=(clusters.size()*xysum-xsum*ysum)/(clusters.size()*x2sum-xsum*xsum); //calculate slope
    b=(x2sum*ysum-xsum*xysum)/(x2sum*clusters.size()-xsum*xsum);           //calculate intercept
    return;
  }


  /**
   * r1 is radius of sPHENIX layer
   * r2, x2 and y2 are parameters of circle fitted to TPC clusters
   * the solutions are xplus, xminus, yplus, yminus

   * The intersection of two circles occurs when
   * (x-x1)^2 + (y-y1)^2 = r1^2,  / (x-x2)^2 + (y-y2)^2 = r2^2
   * Here we assume that circle 1 is an sPHENIX layer centered on x1=y1=0, and circle 2 is arbitrary
   * x^2 +y^2 = r1^2,   (x-x2)^2 + (y-y2)^2 = r2^2
   * expand the equations and subtract to eliminate the x^2 and y^2 terms, gives the radical line connecting the intersection points
   * iy = - (2*x2*x - D) / 2*y2,
   * then substitute for y in equation of circle 1
   */
  bool circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus)
  {

    const double D = square(r1) - square(r2) + square(x2) + square(y2);
    const double a = 1.0 + square(x2/y2);
    const double b = - D * x2/square(y2);
    const double c = square(D/(2.*y2)) - square(r1);
    const double delta = square(b)-4.*a*c;
    if( delta < 0 ) return false;

    const double sqdelta = std::sqrt( delta );

    xplus = (-b + sqdelta ) / (2. * a);
    xminus = (-b - sqdelta ) / (2. * a);

    // both values of x are valid
    // but for each of those values, there are two possible y values on circle 1
    // but only one of those falls on the radical line:

    yplus  = -(2*x2*xplus - D) / (2.*y2);
    yminus = -(2*x2*xminus - D) / (2.*y2);
    return true;

  }

  /// calculate intersection from circle to line, in 2d. return true on success
  /**
   * circle is defined as (x-xc)**2 + (y-yc)**2 = r**2
   * line is defined as nx(x-x0) + ny(y-y0) = 0
   * to solve we substitute y by y0 - nx/ny*(x-x0) in the circle equation and solve the 2nd order polynom
   * there is the extra complication that ny can be 0 (vertical line) to prevent this, we multiply all terms of the polynom by ny**2
   * and account for this special case when calculating x from y
   */
  bool circle_line_intersection(
    double r, double xc, double yc,
    double x0, double y0, double nx, double ny,
    double &xplus, double &yplus, double &xminus, double &yminus)
  {
    if( ny == 0 )
    {
      // vertical lines are defined by ny=0 and x = x0
      xplus = xminus = x0;

      // calculate y accordingly
      const double delta = square(r) - square(x0-xc);
      if( delta < 0 ) return false;

      const double sqdelta = std::sqrt( delta );
      yplus = yc + sqdelta;
      yminus = yc - sqdelta;

    } else {

      const double a = square(nx) + square(ny);
      const double b = -2.*( square(ny)*xc + square(nx)*x0 + nx*ny*(y0-yc) );
      const double c = square(ny)*(square(xc)-square(r)) + square(ny*(y0-yc)+nx*x0);
      const double delta = square(b) - 4.*a*c;
      if( delta < 0 ) return false;

      const double sqdelta = std::sqrt( delta );
      xplus = (-b + sqdelta)/(2.*a);
      xminus = (-b - sqdelta)/(2.*a);

      yplus = y0 -(nx/ny)*(xplus-x0);
      yminus = y0 - (nx/ny)*(xminus-x0);

    }

    return true;

  }

  // streamer of TVector3
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << "," << vector.y() << "," << vector.z() << ")";
    return out;
  }

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

  fdrphi = new TF1("fdrphi", "[0] + [1]*std::abs(x)");
  fdrphi->SetParameter(0, _par0 *_collision_rate / _reference_collision_rate);
  fdrphi->SetParameter(1, _par1);

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

    auto tracklet_tpc = phtrk_iter->second;

    if (Verbosity() >= 1)
    {
      std::cout << std::endl
        << __LINE__
        << ": Processing seed itrack: " << phtrk_iter->first
        << ": nhits: " << tracklet_tpc-> size_cluster_keys()
        << ": Total tracks: " << _track_map->size()
        << ": phi: " << tracklet_tpc->get_phi()
        << std::endl;
    }

    // Get the outermost TPC clusters for this tracklet
    std::map<unsigned int, TrkrCluster*> outer_clusters;
    std::vector<TrkrCluster*> clusters;

    for (SvtxTrack::ConstClusterKeyIter key_iter = tracklet_tpc->begin_cluster_keys(); key_iter != tracklet_tpc->end_cluster_keys(); ++key_iter)
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

    // loop over micromegas layer
    for(unsigned int imm = 0; imm < _n_mm_layers; ++imm)
    {

      // get micromegas geometry object
      const unsigned int layer = _min_mm_layer + imm;
      const auto layergeom = static_cast<CylinderGeomMicromegas*>(_geomContainerMicromegas->GetLayerGeom(layer));
      const auto layer_radius = layergeom->get_radius();

      // method to find where fitted circle intersects this layer
      double xplus = 0;
      double xminus = 0;
      double yplus = 0;
      double yminus = 0;

      // finds the intersection of the fitted circle with the micromegas layer
      if( !circle_circle_intersection(	layer_radius, R, X0, Y0, xplus, yplus, xminus, yminus) )
      {
        if(Verbosity() > 10)
        {
          std::cout << PHWHERE << " circle/circle intersection calculation failed, skip this case" << std::endl;
          std::cout << PHWHERE << " mm_radius " << layer_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
        }

        continue;
      }

      // we can figure out which solution is correct based on the last cluster position in the TPC
      const double last_clus_phi = std::atan2(clusters.back()->getY(), clusters.back()->getX());
      double phi_plus = std::atan2(yplus, xplus);
      double phi_minus = std::atan2(yminus, xminus);

      // calculate z
      double r = layer_radius;
      double z = B + A*r;

      // select the angle that is the closest to last cluster
      // store phi, apply coarse space charge corrections in calibration mode
      double phi = std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus) ? phi_plus:phi_minus;
      if(_sc_calib_mode) phi -= fdrphi->Eval(z)/r;

      // create cylinder intersection point in world coordinates
      const TVector3 world_intersection_cylindrical( r*std::cos(phi), r*std::sin(phi), z );

      // find matching tile
      int tileid = layergeom->find_tile_cylindrical( world_intersection_cylindrical );
      if( tileid < 0 ) continue;

      // get tile coordinates
      const auto tile_center_phi = layergeom->get_center_phi( tileid );
      const double x0 = r*std::cos( tile_center_phi );
      const double y0 = r*std::sin( tile_center_phi );

      const double nx = x0;
      const double ny = y0;

      // calculate intersection to tile
      if( !circle_line_intersection( R, X0, Y0, x0, y0, nx, ny, xplus, yplus, xminus, yminus) )
      {
        if(Verbosity() > 10)
        { std::cout << PHWHERE << "circle_line_intersection - failed" << std::endl; }
        continue;
      }

      // select again angle closest to last cluster
      phi_plus = std::atan2(yplus, xplus);
      phi_minus = std::atan2(yminus, xminus);
      const bool is_plus = (std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus));

      // calculate x, y and z
      const double x = (is_plus ? xplus:xminus);
      const double y = (is_plus ? yplus:yminus);
      r = get_r( x, y );
      z = B + A*r;

      /*
       * create planar intersection point in world coordinates
       * this is the position to be compared to the clusters
       */
      const TVector3 world_intersection_planar( x, y, z );

      // convert to tile local reference frame, apply SC correction
      TVector3 local_intersection_planar = layergeom->get_local_from_world_coords( tileid, world_intersection_planar );
      if(_sc_calib_mode)
      {
        /*
         * apply SC correction to the local intersection,
         * to make sure that the corrected intersection is still in the micromegas plane
         * in local tile coordinates, the rphi direction, to which the correction is applied, corresponds to the x direction
         */
        local_intersection_planar.SetX( local_intersection_planar.x() - fdrphi->Eval(z) );
      }

      // store segmentation type
      const auto segmentation_type = layergeom->get_segmentation_type();

      // generate tilesetid and get corresponding clusters
      const auto tilesetid = MicromegasDefs::genHitSetKey(layer, segmentation_type, tileid);
      const auto mm_clusrange = _cluster_map->getClusters(tilesetid);

      // convert to tile local coordinate and compare
      for(auto clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
      {

        // store cluster and key
        const auto& [key, cluster] = *clusiter;

        const TVector3 world_cluster( cluster->getX(), cluster->getY(), cluster->getZ() );
        const TVector3 local_cluster = layergeom->get_local_from_world_coords( tileid, world_cluster );

        // compute residuals and store
        /* in local tile coordinate, x is along rphi, and z is along z) */
        const double drphi = local_intersection_planar.x() - local_cluster.x();
        const double dz = local_intersection_planar.z() - local_cluster.z();

        // compare to cuts and add to track if matching
        if( std::abs(drphi) < _rphi_search_win[imm] && std::abs(dz) < _z_search_win[imm] )
        {
          tracklet_tpc->insert_cluster_key(key);
          _assoc_container->SetClusterTrackAssoc(key, tracklet_tpc->get_id());

          // prints out a line that can be grep-ed from the output file to feed to a display macro
          if( _test_windows )
          {
            // cluster rphi and z
            const double mm_clus_rphi = get_r( cluster->getX(), cluster->getY() ) * std::atan2( cluster->getY(),  cluster->getX() );
            const double mm_clus_z = cluster->getZ();

            // projection phi and z, without correction
            const double rphi_proj = get_r( world_intersection_planar.x(), world_intersection_planar.y() ) * std::atan2( world_intersection_planar.y(), world_intersection_planar.x() );
            const double z_proj = world_intersection_planar.z();

            /*
             * Note: drphi and dz might not match the difference of the rphi and z quoted values. This is because
             * 1/ drphi and dz are actually calculated in Tile's local reference frame, not in world coordinates
             * 2/ drphi also includes SC distortion correction, which the world coordinates don't
             */
            std::cout
              << "  Try_mms: " << (int) layer
              << " drphi " << drphi
              << " dz " << dz
              << " mm_clus_rphi " << mm_clus_rphi << " mm_clus_z " << mm_clus_z
              << " rphi_proj " <<  rphi_proj << " z_proj " << z_proj
              << std::endl;
          }
        }

      } // end loop over clusters

    } // end loop over Micromegas layers

    if(Verbosity() > 3)
    { tracklet_tpc->identify(); }

  }

  if(Verbosity() > 0)
  { std::cout << " Final track map size " << _track_map->size() << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________________________
int PHMicromegasTpcTrackMatching::End()
{ return Fun4AllReturnCodes::EVENT_OK; }

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

  return Fun4AllReturnCodes::EVENT_OK;
}

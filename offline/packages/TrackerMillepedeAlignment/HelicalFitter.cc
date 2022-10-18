#include "HelicalFitter.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase/TrkrClusterContainer.h>   
#include <trackbase/TrkrClusterCrossingAssoc.h>   
#include <trackbase/TrackFitUtils.h>

#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

using namespace std;

//____________________________________________________________________________..
HelicalFitter::HelicalFitter(const std::string &name):
  SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//____________________________________________________________________________..
HelicalFitter::~HelicalFitter()
{

}

//____________________________________________________________________________..
int HelicalFitter::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();

   int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//_____________________________________________________________________
void HelicalFitter::SetDefaultParameters()
{


  return;
}

//____________________________________________________________________________..
int HelicalFitter::process_event(PHCompositeNode*)
{
  // _track_map_tpc contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

  if(Verbosity() > 0)
    cout << PHWHERE 
	 << " TPC track map size " << _track_map_tpc->size() 
	 << " Silicon track map size "  << _track_map_silicon->size() 
	 << endl;

  if(_track_map_silicon->size() == 0 && _track_map_tpc->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;

  // Decide whether we want to make a helical fit for silicon or TPC
  unsigned int maxtracks = 0; 
  if(fitsilicon)  { maxtracks =  _track_map_silicon->size(); }
  else if(fittpc) { maxtracks =  _track_map_tpc->size();  }

  for(unsigned int trackid = 0; trackid < maxtracks; ++trackid)
    {
      TrackSeed* _tracklet = nullptr;
      if(fitsilicon) {  _tracklet = _track_map_silicon->get(trackid); }
      else if(fittpc) {  _tracklet = _track_map_tpc->get(trackid);	 }

      if(!_tracklet)
	{
	  trackid++;
	  continue;
	}

      // The seed fit parameters are for uncorrected clusters, so we fit corrected positions

      // loop over all clusters
      std::vector<Acts::Vector3> global_vec;
      std::vector<TrkrDefs::cluskey> cluskey_vec;

      for (auto clusIter = _tracklet->begin_cluster_keys();
	   clusIter != _tracklet->end_cluster_keys();
	   ++clusIter)
	{
	  auto key = *clusIter;
	  auto cluster = _cluster_map->findCluster(key);
	  if(!cluster)
	    {
	      std::cout << "Failed to get cluster with key " << key << " for track " << trackid << std::endl;
	      continue;
	    }	  
	  
	  //	  auto subsurfkey = cluster->getSubSurfKey();
	  
	  /// Make a safety check for clusters that couldn't be attached to a surface
	  auto surf = _tGeometry->maps().getSurface(key, cluster);
	  if(!surf)
	    { continue; }
	  
	  Acts::Vector3 global  = _tGeometry->getGlobalPosition(key, cluster);	  

	  const unsigned int trkrid = TrkrDefs::getTrkrId(key);	  
	  if(trkrid ==  TrkrDefs::tpcId)
	    {	  
	      // make all corrections to global position of TPC cluster
	      const unsigned int side = TpcDefs::getSide(key);
	      int crossing = 0;
	      float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
	      global[2] = z;
	      
	      // apply distortion corrections
	      if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
	      if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
	      if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
	    }
	  
	  // add the global positions to a vector to give to the helical fitter
	  global_vec.push_back(global);
	  cluskey_vec.push_back(key);
      	  
	} // end loop over clusters for this track

      // make the helical fit using TrackFitUtils
      std::tuple<double, double, double> circle_fit_pars = TrackFitUtils::circle_fit_by_taubin(global_vec);
      std::tuple<double,double> line_fit_pars = TrackFitUtils::line_fit(global_vec);

      // capture fit pars, put on node tree
      float radius = std::get<0>(circle_fit_pars);
      float x0 = std::get<1>(circle_fit_pars);
      float y0 = std::get<2>(circle_fit_pars);
      float zslope = std::get<0>(line_fit_pars);
      float z0 = std::get<1>(line_fit_pars);
      Acts::Vector3 helix_center(x0,y0,z0);        //  center of circle in xy plane, and starting z position

      std::cout << " Track " << trackid << " helix center and start z " << x0 << " " << y0 << " " << z0 
		<< " radius " << radius << " zslope " << zslope << std::endl;

      // get the residuals and derivatives for all clusters
      std::vector<Acts::Vector3> residual_vec;
      std::vector<float> derivative_vec;
      for( const auto& global : global_vec )
	{
	  // PCA of helix to cluster global position
	  Acts::Vector3 pca = get_helix_pca(radius, zslope, helix_center, global);
	    
	  // capture residuals
	  residual_vec.push_back(global - pca);

	  std::cout << "    cluster position " << global(0) << " " << global(1) << " " << global(2) 
		    << " pca " << pca(0) << " " << pca(1) << " " << pca(2)  << std::endl;

	  // need the (non-zero) derivatives with respect to the track parameters
	  // parameters are radius, circle (x0,y0), z0 and zslope

	  // z = z0 + zslope * radius 
	  float z_zslope = radius;
	  float z_z0 = 1.0;
	  derivative_vec.push_back(z_zslope);
	  derivative_vec.push_back(z_z0);
	  //float z_x = 0.0;  // could skip this
	  //float z_y = 0.0;  // could skip this
	  //derivative_vec.push_back(z_x);
	  //derivative_vec.push_back(z_y);
	  //std::cout << "    z derivatives: z_zslope " << z_zslope << " z_z0 " << z_z0 << " z_x " << z_x << " z_y " << z_y << std::endl;
	  std::cout << "    z derivatives: z_zslope " << z_zslope << " z_z0 " << z_z0 << std::endl;

	  // dx/dradius = radius / (R^2 - (y-y0)^2)^1/2  , sign of (x-x0)/fabs(x-x0)
	  float x_x0 = 1.0;
	  float x_radius = ( pca(0)-helix_center(0) / fabs((pca(0)-helix_center(0)) ) ) * radius / sqrt(radius*radius - (pca(1) - helix_center(1) * (pca(1) - helix_center(1) ) ) );
	  derivative_vec.push_back(x_x0);
	  derivative_vec.push_back(x_radius);
	  //float x_y = 0.0;  // could skip this
	  //float x_z = 0.0;  // could skip this
	  //derivative_vec.push_back(x_y);
	  //derivative_vec.push_back(x_z);
	  //std::cout << "    x derivatives: x_x0 " << x_x0 << " x_radius " << x_radius << " x_y " << x_y << " x_z " << x_z << std::endl;
	  std::cout << "    x derivatives: x_x0 " << x_x0 << " x_radius " << x_radius << std::endl;

	  // dy/dradius = sign * radius / (R^2 - (x-x0)^2)^1/2 , sign of (y-y0)/fabs(y-y0)
	  float y_y0 = 1.0;
	  float y_radius =  ( pca(1)-helix_center(1) / fabs((pca(1)-helix_center(1)) ) ) * radius / sqrt(radius*radius - (pca(0) - helix_center(0) * (pca(0) - helix_center(0) ) ) );
	  derivative_vec.push_back(y_y0);
	  derivative_vec.push_back(y_radius);
	  //float y_x = 0.0;  // could skip this
	  //float y_z = 0.0;  // could skip this
	  //derivative_vec.push_back(y_x);
	  //derivative_vec.push_back(y_z);
	  //std::cout << "    y derivatives: y_y0 " << y_y0 << " y_radius " << y_radius << " y_x " << y_x << " y_z " << y_z << std::endl;
	  std::cout << "    y derivatives: y_y0 " << y_y0 << " y_radius " << y_radius << std::endl;

	  // now we need the derivatives wrt the alignment parameters
	  // copy methods over from MakeMilleFiles? Or write node tree info to pick up by MakeMilleFiles?



	  // add everything to the Mille file for this cluster


	}

    }  // end loop over tracks
  



  return Fun4AllReturnCodes::EVENT_OK;

 }
  
int HelicalFitter::End(PHCompositeNode* )
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  HelicalFitter::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_tpc = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map_tpc)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
    {
      std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
} 

Acts::Vector3 HelicalFitter::get_helix_pca(float radius, float zslope, Acts::Vector3 helix_center, Acts::Vector3 global)
{
  // no analytic solution for the coordinates of the closest approach of a helix to a point
  // Instead, we get the PCA in x and y to the circle, and the PCA in z to the z vs R line at the R of the PCA 
  
  Acts::Vector2 pca_circle = get_circle_point_pca(radius, helix_center, global);

  // The radius of the PCA determines the z position:
  float pca_circle_radius = pca_circle.norm();
  float pca_z = pca_circle_radius * zslope + helix_center(2);
  Acts::Vector3 pca(pca_circle(0), pca_circle(1), pca_z);

  // now we want a second point on the helix so we can get a local straight line approximation to the track
  // project the circle PCA vector an additional small amount and find the helix PCA to that point 
  float projection = 0.25;  // cm
  Acts::Vector3 second_point = pca + projection * pca/pca.norm();
  Acts::Vector2 second_point_pca_circle = get_circle_point_pca(radius, helix_center, second_point);
  float second_point_pca_z = pca_circle_radius * zslope + helix_center(2);
  Acts::Vector3 second_point_pca(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) /  (second_point_pca - pca).norm();

 // get the PCA of the cluster to that line
  Acts::Vector3 final_pca = getPCALinePoint(global, tangent, pca);

  return final_pca;
}

Acts::Vector3 HelicalFitter::getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)   

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ( (global - posref).dot(tangent) ) * tangent;

  return pca;
}

Acts::Vector2 HelicalFitter::get_circle_point_pca(float radius, Acts::Vector3 center, Acts::Vector3 global)
{
  // get the PCA of a cluster (x,y) position to a circle
  // draw a line from the origin of the circle to the point
  // the intersection of the line with the circle is at the distance radius from the origin along that line 

  Acts::Vector2 origin(center(0), center(1));
  Acts::Vector2 point(global(0), global(1));

  Acts::Vector2 pca = origin + radius * (point - origin) / (point - origin).norm();

  return pca;
}


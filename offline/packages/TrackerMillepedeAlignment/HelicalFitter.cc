#include "HelicalFitter.h"

#include "Mille.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterContainer.h>   
#include <trackbase/TrkrClusterCrossingAssoc.h>   

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
  , _mille(nullptr)
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

  // Instantiate Mille and open output data file
  if(test_output)
    {
      _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  }
 else
   {
     _mille = new Mille(data_outfilename.c_str()); 
   }

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
  steering_file.close();

  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are silicon " << si_grp << " tpc " << tpc_grp << " mms " << mms_grp << std::endl; 

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
	 << " TPC seed map size " << _track_map_tpc->size() 
	 << " Silicon seed map size "  << _track_map_silicon->size() 
	 << endl;

  if(_track_map_silicon->size() == 0 && _track_map_tpc->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;

  // Decide whether we want to make a helical fit for silicon or TPC
  unsigned int maxtracks = 0; 
  if(fittpc) { maxtracks =  _track_map_tpc->size();  }
  if(fitsilicon)  { maxtracks =  _track_map_silicon->size(); }
  for(unsigned int trackid = 0; trackid < maxtracks; ++trackid)
    {
      TrackSeed* tracklet = nullptr;
      if(fitsilicon) {  tracklet = _track_map_silicon->get(trackid); }
      else if(fittpc) {  tracklet = _track_map_tpc->get(trackid);	 }
      if(!tracklet) { trackid++;  continue; }

      std::vector<Acts::Vector3> global_vec;
      std::vector<TrkrDefs::cluskey> cluskey_vec;


      ////getTrackletClusters(tracklet, global_vec, cluskey_vec);   // store cluster corrected global positions in a vector 

      // Get a vector of cluster keys from the tracklet  
      getTrackletClusterList(tracklet, cluskey_vec);
      // store cluster global positions in a vector
      TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   

      correctTpcGlobalPositions( global_vec, cluskey_vec);

      ////std::vector<float> fitpars =  fitClusters(global_vec, cluskey_vec);       // do helical fit
      std::vector<float> fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit`

      if(fitpars.size() == 0) continue;  // discard this track, not enough clusters to fit

      if(Verbosity() > 1)  
	{ std::cout << " Track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
		 << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }

      // if a full track is requested, get the silicon clusters too and refit
      if(fittpc && fitfulltrack)
	{
	  // this associates silicon clusters and adds them to the vectors
	  //unsigned int nsilicon = addSiliconClusters(fitpars, global_vec, cluskey_vec);
	  unsigned int nsilicon = TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
	  if(nsilicon < 3) continue;  // discard this TPC seed, did not get a good match to silicon

	  // fit the full track now
	  fitpars.clear();
	  ////fitpars =  fitClusters(global_vec, cluskey_vec);       // do helical fit
	  fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
	  if(fitpars.size() == 0) continue;  // discard this track, fit failed

	  if(Verbosity() > 1)  
	    { std::cout << " Full track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
					   << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }
	} 

      // get the residuals and derivatives for all clusters
      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{
	  auto global = global_vec[ivec];
	  auto cluskey = cluskey_vec[ivec];
	  auto cluster = _cluster_map->findCluster(cluskey);
	  if(!cluster) { continue;}

	  // What we need now is to find the point on the surface at which the helix would intersect
	  // If we have that point, we can transform the fit back to local coords
	  // we have fitpars for the helix and the cluster key, from which we get the surface

	  Surface surf = _tGeometry->maps().getSurface(cluskey, cluster);
	  Acts::Vector3 fitpoint = get_helix_surface_intersection(surf, fitpars, global);
	  std::cout << " fitpoint " << fitpoint(0) << "  " << fitpoint(1) << "  " << fitpoint(2) << std::endl;

	  // fitpoint is the point where the helical fit intersects the plane of the surface
	  // this is what we need to get the residuals
       	    
	  // Now transform the helix fitpoint to local coordinates to compare with cluster local coordinates
	  Acts::Vector3 fitpoint_local = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (fitpoint *  Acts::UnitConstants::cm);
	  fitpoint_local /= Acts::UnitConstants::cm;

	  Acts::Vector2 residual(cluster->getLocalX() - fitpoint_local(0), cluster->getLocalY() - fitpoint_local(1)); 
  
	  if(Verbosity() > 1) {std::cout << "    cluster position " << global(0) << " " << global(1) << " " << global(2) << std::endl
					 << " fitpoint " << fitpoint(0) << " " << fitpoint(1) << " " << fitpoint(2) << std::endl
					 << " fitpoint_local " << fitpoint_local(0) << " " << fitpoint_local(1) << " " << fitpoint_local(2) << std::endl  
					 << " cluster local x " << cluster->getLocalX() << " cluster local y " << cluster->getLocalY() << std::endl
					 << " cluster local residual x " << residual(0) << " cluster local residual y " <<residual(1) << std::endl;}
	  unsigned int layer = TrkrDefs::getLayer(cluskey_vec[ivec]);	  
	  float phi =  atan2(global(1), global(0));
	  std::cout << "Local residuals: layer " << layer << " phi " << phi << " dx " << residual(0) << " dy " << residual(1) << std::endl;

	  if(Verbosity() > 2)
	    {
	      if(layer < 2)
		{
		  std::cout << " Global residuals: layer " << layer << " phi " << phi
			    << " dx " << global(0) - fitpoint(0) 
			    << " dy " << global(1) - fitpoint(1) 
			    << std::endl;
		  /*
		  Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitpars, global);
		  std::cout << " layer " << layer << " phi " << phi
			    << " dx " << global(0) - pca(0) 
			    << " dy " << global(1) - pca(1) 
			    << std::endl;
		  */
		  /*		  
		  // check
		  float phirel = atan2( (global(0) - fitpars[1]), (global(1) - fitpars[2]));
		  float pdist = sqrt( pow(global(0) - fitpars[1], 2) + pow(global(1) - fitpars[2], 2));
		  float pcadist = sqrt( pow(fitpoint(0) - fitpars[1], 2) + pow(fitpoint(1) - fitpars[2], 2));
		  std::cout << " helix fitpoint:layer " << layer << " phi " << phi << " phirel " << phirel << " pdist " << pdist << " pcadist " << pcadist 
			    << " R " << fitpars[0] << std::endl; 
		  */
		}
	    }

	  // need standard deviation of measurements
	  Acts::Vector3 clus_sigma = getClusterError(cluster, cluskey, global);
	  if(isnan(clus_sigma(0)) || isnan(clus_sigma(1)) || isnan(clus_sigma(2)))  { continue; }

	  int glbl_label[NGL];
	  getGlobalLabels(surf, glbl_label);  // these depend on the sensor grouping

	  float lcl_derivative[NLC];
	  float glbl_derivative[NGL];
	  // The angleDerivs dimensions are [alpha/beta/gamma](x/y/z)
	  std::vector<Acts::Vector3> angleDerivs = getDerivativesAlignmentAngles(global, cluskey, cluster); 
	  std::vector<Acts::Vector3> translDerivs = getDerivativesAlignmentTranslations(global, cluskey, cluster);

	  // Add the measurement separately for each coordinate direction to Mille
	  // set the derivatives non-zero only for parameters we want to be optimized
	  getLocalDerivativesX(surf, fitpoint, fitpoint_local, fitpars, lcl_derivative);
	  getGlobalDerivativesX(angleDerivs, translDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3)
	    { std::cout << "layer " << layer << " X buffers:" << std::endl; printBuffers(0, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if( !isnan(residual(0)) && clus_sigma(0) < 1.0)  // discards crazy clusters
	    { _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(0), _error_inflation*clus_sigma(0));}

	  /*
	  getLocalDerivativesY(surf, fitpoint, fitpoint_local, fitpars, lcl_derivative);
	  getGlobalDerivativesY(angleDerivs, translDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3) 
	    { std::cout  << "layer " << layer << " Y buffers:" << std::endl; printBuffers(1, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if( !isnan(residual(1)) && clus_sigma(1) < 1.0)  // discards crazy clusters
	    {_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(1), _error_inflation*clus_sigma(1));}
	  */
	  
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
	  getLocalDerivativesZ(fitpoint, fitpars, lcl_derivative);
	  getGlobalDerivativesZ(angleDerivs, translDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3) 
	    { std::cout  << "layer " << layer << " Z buffers:" << std::endl; printBuffers(2, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if(!isnan(residual(1)) && clus_sigma(2) < 1.0 && trkrid != TrkrDefs::inttId)
	    {_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(1), _error_inflation*clus_sigma(2));}

	}

      // close out this track
      _mille->end();
      
    }  // end loop over tracks
  
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 HelicalFitter::get_helix_surface_intersection(Surface surf, std::vector<float>& fitpars, Acts::Vector3 global)
{
  // we want the point where the helix intersects the plane of the surface

  // get the plane of the surface
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
  Acts::Vector3 sensorNormal    = -surf->normal(_tGeometry->geometry().getGeoContext());
  sensorNormal /= sensorNormal.norm();

  // there are analytic solutions for a line-plane intersection.
  // to use this, need to get the vector tangent to the helix near the measurement and a point on it.
  std::pair<Acts::Vector3, Acts::Vector3> line =  get_helix_tangent(fitpars, global);
  Acts::Vector3 pca = line.first;
  Acts::Vector3 tangent = line.second;

  std::cout << "   pca: "  << pca(0) << "  " << pca(1) << "  " << pca(2) << "  " << std::endl;

  Acts::Vector3 intersection = get_line_plane_intersection(pca, tangent, sensorCenter, sensorNormal);

  return intersection;
}

Acts::Vector3 HelicalFitter::get_line_plane_intersection(Acts::Vector3 PCA, Acts::Vector3 tangent, Acts::Vector3 sensor_center, Acts::Vector3 sensor_normal)
{
  // get the intersection of the line made by PCA and tangent with the plane of the sensor

  // For a point on the line
  // p = PCA + d * tangent;
  // for a point on the plane
  // (p - sensor_center).sensor_normal = 0

  // The solution is:
  float d = (sensor_center - PCA).dot(sensor_normal) / tangent.dot(sensor_normal);
  Acts::Vector3 intersection = PCA + d * tangent;
  std::cout << "   intersection: " << intersection(0) << "  "  << intersection(1) << "  "  << intersection(2) << "  " << std::endl;
  std::cout << "        sensor_center: " << sensor_center(0) << "  " << sensor_center(1) << "  " << sensor_center(2) << "  " << std::endl;
  std::cout << "        sensor_normal: " << sensor_normal(0) << "  " << sensor_normal(1) << "  " << sensor_normal(2) << "  " << std::endl;

  return intersection;
}


std::pair<Acts::Vector3, Acts::Vector3> HelicalFitter::get_helix_tangent(std::vector<float>& fitpars, Acts::Vector3 global)
{
  // no analytic solution for the coordinates of the closest approach of a helix to a point
  // Instead, we get the PCA in x and y to the circle, and the PCA in z to the z vs R line at the R of the PCA 

  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];  
  float zslope = fitpars[3];
  float z0 = fitpars[4];

  Acts::Vector2 pca_circle = get_circle_point_pca(radius, x0, y0, global);

  // The radius of the PCA determines the z position:
  float pca_circle_radius = pca_circle.norm();
  float pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 pca(pca_circle(0), pca_circle(1), pca_z);

  // now we want a second point on the helix so we can get a local straight line approximation to the track
  // project the circle PCA vector an additional small amount and find the helix PCA to that point 
  float projection = 0.25;  // cm
  Acts::Vector3 second_point = pca + projection * pca/pca.norm();
  Acts::Vector2 second_point_pca_circle = get_circle_point_pca(radius, x0, y0, second_point);
  float second_point_pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 second_point_pca(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) /  (second_point_pca - pca).norm();

 // get the PCA of the cluster to that line
  Acts::Vector3 final_pca = getPCALinePoint(global, tangent, pca);

  std::pair<Acts::Vector3, Acts::Vector3> line = std::make_pair(final_pca, tangent);

  return line;
}
  
int HelicalFitter::End(PHCompositeNode* )
{
  // closes output file in destructor
  delete _mille;

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

Acts::Vector3 HelicalFitter::get_helix_pca(std::vector<float>& fitpars, Acts::Vector3 global)
{
  return TrackFitUtils::get_helix_pca(fitpars, global);
}


Acts::Vector3 HelicalFitter::getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)   

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ( (global - posref).dot(tangent) ) * tangent;

  return pca;
}

Acts::Vector2 HelicalFitter::get_circle_point_pca(float radius, float x0, float y0, Acts::Vector3 global)
{
  // get the PCA of a cluster (x,y) position to a circle
  // draw a line from the origin of the circle to the point
  // the intersection of the line with the circle is at the distance radius from the origin along that line 

  Acts::Vector2 origin(x0, y0);
  Acts::Vector2 point(global(0), global(1));

  Acts::Vector2 pca = origin + radius * (point - origin) / (point - origin).norm();

  return pca;
}

std::vector<Acts::Vector3> HelicalFitter::getDerivativesAlignmentTranslations(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  std::vector<Acts::Vector3> derivs_vector;

  // Make a transform that applies small translations in the surface frame
  for(unsigned int itrans = 0; itrans < 3; ++itrans)
    {
      // creates transform that adds a perturbation translation along one axis, uses it to estimate derivative wrt perturbation translation

      unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      Acts::Vector3 derivs(0,0,0);
      Eigen::Vector3d theseTransl(0,0,0);
      theseTransl[itrans] = sensorTransl[itrans];  // set the one we want to be non-zero

      Acts::Vector3 keeper(0,0,0);
      for(int ip = 0; ip < 2; ++ip)
	{
	  if(ip == 1) { theseTransl[itrans] *= -1; } // test both sides of zero

	  if(Verbosity() > 1)
	    { std::cout << "     trkrId " << trkrId << " layer " << layer << " cluster_key " << cluster_key 
			<< " sensorTransl " << theseTransl[0] << "  " << theseTransl[1] << "  " << theseTransl[2] << std::endl; }

	  Acts::Transform3 perturbationTranslation = makePerturbationTranslation(theseTransl);
	  
	  // transform the cluster local position to global coords with this additional translation added
	  auto x = cluster->getLocalX() * 10;  // mm 
	  auto y = cluster->getLocalY() * 10;	  
	  if(trkrId == TrkrDefs::tpcId) { y = convertTimeToZ(cluster_key, cluster); }
	  
	  Eigen::Vector3d clusterLocalPosition (x,0,y);  // follows our convention for local coords
	  Eigen::Vector3d finalCoords = perturbationTranslation*clusterLocalPosition;  // result in mm
	  finalCoords /= 10.0; // convert mm back to cm

	  // have to add corrections for TPC clusters after transformation to global
	  //	  if(trkrId == TrkrDefs::tpcId) {  makeTpcGlobalCorrections(cluster_key, crossing, global); }

	  // note that x and y cancel out here	  	  
	  if(ip == 0)
	    {
	      keeper(0) = (finalCoords(0) - x);
	      keeper(1) = 0;
	      keeper(2) = (finalCoords(2) - y);
	    }
	  else
	    {
	      keeper(0) -= (finalCoords(0) - x);
	      keeper(1) -= 0;
	      keeper(2) -= (finalCoords(2) - y);
	    }

	  if(Verbosity() > 5)
	    {
	      std::cout << "        AlignmentTranslationsDerivs: finalCoords(0) " << finalCoords(0) << " global(0) " << global(0) << " finalCoords(1) " 
			<< finalCoords(1) << " global(1) " << global(1) << " finalCoords(2) " << finalCoords(2) 
			<< " global(2) " << global(2) << std::endl;
	      std::cout  << "        keeper now:  keeper(0) " << keeper(0) << " keeper(1) " << keeper(1) << " keeper(2) " 
			 << keeper(2) << std::endl;
	    }
	}

      // derivs vector contains:
      //   (dx/dx,     dy/dx,     dz/dx)    (for itrans = 0)
      //   (dx/dy,      dy/dy,     dz/dy)   (for itrans = 1) 
      //   (dx/dz,      dy/dz,     dz/dz)   (for itrans = 2)
 
     // Average the changes to get the estimate of the derivative      
      derivs(0) = keeper(0) / (2.0 * 0.1 * fabs(theseTransl[itrans]));  // convert theseTransl to cm 
      if( isnan(derivs(0)) ) { derivs(0) = 0; }
      derivs(1) = keeper(1) / (2.0 * 0.1 * fabs(theseTransl[itrans]));
      if( isnan(derivs(1)) ) { derivs(1) = 0; }
      derivs(2) = keeper(2) / (2.0 * 0.1 * fabs(theseTransl[itrans]));
      if( isnan(derivs(2)) ) { derivs(2) = 0; }
      derivs_vector.push_back(derivs);

      if(Verbosity() > 1) { std::cout << "        AlignmentTranslationsDerivs: itrans " << itrans << " derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl; }
    }
  return derivs_vector;
}

std::vector<Acts::Vector3> HelicalFitter::getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  // Wev want the effect of  a small rotation around the relevant axis in the local frame on the local coords

  std::vector<Acts::Vector3> derivs_vector;

  for(unsigned int iangle = 0; iangle < 3; ++iangle)
    {
      // creates transform that adds a perturbation rotation around one axis, uses it to estimate derivative wrt perturbation rotation

      unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      Acts::Vector3 derivs(0,0,0);
      Eigen::Vector3d theseAngles(0,0,0);
      theseAngles[iangle] = sensorAngles[iangle];  // set the one we want to be non-zero

      Acts::Vector3 keeper(0,0,0);
      for(int ip = 0; ip < 2; ++ip)
	{
	  if(ip == 1) { theseAngles[iangle] *= -1; } // test both sides of zero

	  if(Verbosity() > 1)
	    { std::cout << "     trkrId " << trkrId << " layer " << layer << " cluster_key " << cluster_key 
			<< " sensorAngles " << theseAngles[0] << "  " << theseAngles[1] << "  " << theseAngles[2] << std::endl; }

	  Acts::Transform3 perturbationTransformation = makePerturbationTransformation(theseAngles);

	  // transform the cluster local position using this additional rotation
	  auto x = cluster->getLocalX() * 10;  // mm 
	  auto y = cluster->getLocalY() * 10;	  
	  if(trkrId == TrkrDefs::tpcId) { y = convertTimeToZ(cluster_key, cluster); }
	  
	  Eigen::Vector3d clusterLocalPosition (x,0,y);  // our convention, applies until Acts global rotation occurs
	  Eigen::Vector3d finalLocalPosition = perturbationTransformation*clusterLocalPosition;  // result in mm
	  finalLocalPosition /= 10.0; // convert mm back to cm

	  // have to add corrections for TPC clusters after transformation to global
	  // The helical fit is to corrected data, so if we transform back to local, can 
	  // we compare with the cluster local? What is needed here for the TPC?

	  // note that x and y cancel out here	  	  
	  if(ip == 0)
	    {
	      keeper(0) = (finalLocalPosition(0) - x);
	      keeper(2) = (finalLocalPosition(2) - y);
	    }
	  else
	    {
	      keeper(0) -= (finalLocalPosition(0) - x);
	      keeper(2) -= (finalLocalPosition(2) - y);
	    }

	  if(Verbosity() > 5)
	    {
	      std::cout << "        AlignmentAngleDerivs: finalLocalPosition(0) " << finalLocalPosition(0) << " global(0) " << global(0) << " finalLocalPosition(1) " << finalLocalPosition(1) 
			<< " global(1) " << global(1) << " finalLocalPosition(2) " << finalLocalPosition(2) << " global(2) " << global(2) << std::endl;
	      std::cout  << "        keeper now:  keeper(0) " << keeper(0) << " keeper(1) " << keeper(1) << " keeper(2) " << keeper(2) << std::endl;
	    }
	}

      // derivs vector contains:
      //   (dx/dalpha,   dz/dalpha)     (for iangle = 0)
      //   (dx/dbeta,      dz/dbeta)        (for iangle = 1) 
      //   (dx/dgamma, dz/dgamma) (for iangle = 2)
 
     // Average the changes to get the estimate of the derivative      
      derivs(0) = keeper(0) / (2.0 * fabs(theseAngles[iangle]));
      if( isnan(derivs(0)) ) { derivs(0) = 0; }
      derivs(1) = 0.0;  // this would be the unused y axis in the local coord frame
      derivs(2) = keeper(2) / (2.0 * fabs(theseAngles[iangle]));
      if( isnan(derivs(2)) ) { derivs(2) = 0; }
      derivs_vector.push_back(derivs);

      if(Verbosity() > 1) { std::cout << "        AlignmentAngleDerivs: iangle " << iangle << " derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl; }
    }
  
  return derivs_vector;
}

  Acts::Transform3 HelicalFitter::makePerturbationTranslation(Acts::Vector3 translations)
  {
    // combine unit rotation and perturbation translations in the local frame into an affine matrix
    Acts::Transform3 perturbationTransformation;

    Eigen::AngleAxisd alpha(0.0, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd beta(0.0, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd gamma(0.0, Eigen::Vector3d::UnitZ());
    Eigen::Quaternion<double> q       = gamma*beta*alpha;
    Eigen::Matrix3d nullRotation = q.matrix();
    perturbationTransformation.linear() = nullRotation;

    perturbationTransformation.translation() = translations;
    
    return perturbationTransformation;    
  }

  Acts::Transform3 HelicalFitter::makePerturbationTransformation(Acts::Vector3 angles)
  {
    // Note: Here beta is applied to the z axis and gamma is applied to the y axis because the geocontext transform 
    // will flip those axes when transforming to global coordinates
    Eigen::AngleAxisd alpha(angles(0), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd beta(angles(1), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd gamma(angles(2), Eigen::Vector3d::UnitZ());
    Eigen::Quaternion<double> q       = gamma*beta*alpha;
    Eigen::Matrix3d perturbationRotation = q.matrix();
    
    // combine rotation and translation into an affine matrix
    Eigen::Vector3d nullTranslation(0,0,0);
    Acts::Transform3 perturbationTransformation;
    perturbationTransformation.linear() = perturbationRotation;
    perturbationTransformation.translation() = nullTranslation;
    
    return perturbationTransformation;    
  }

float HelicalFitter::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = _tGeometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89; // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;   // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if(side == 0) zloc = -zloc;
  float z = zloc * 10.0;
 
  return z; 
}

void HelicalFitter::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;
  
  // apply distortion corrections
  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
}

void HelicalFitter::getGlobalLabels(Surface surf, int glbl_label[])
{
  // identify the global alignment parameters for this surface
  Acts::GeometryIdentifier id = surf->geometryId();
  int label_base = getLabelBase(id);   // This value depends on how the surfaces are grouped	  
  for(int i=0;i<NGL;++i) 
    {
      glbl_label[i] = label_base + i;
      if(Verbosity() > 1)
	{ std::cout << "    glbl " << i << " label " << glbl_label[i] << " "; }
    }
  if(Verbosity() > 1) { std::cout << std::endl; }
}

int HelicalFitter::getTpcRegion(int layer)
{
  int region = 0;
  if(layer > 23 && layer < 39)
    region = 1;
  if(layer > 38 && layer < 55)
    region = 2;

  return region;  
}

int HelicalFitter::getLabelBase(Acts::GeometryIdentifier id)
{
  unsigned int volume = id.volume(); 
  unsigned int acts_layer = id.layer();
  unsigned int layer = base_layer_map.find(volume)->second + acts_layer / 2 -1;
  unsigned int sensor = id.sensitive() - 1;  // Acts starts at 1

  int label_base = 1;  // Mille wants to start at 1

  // decide what level of grouping we want
  if(layer < 7)
    {
      if(si_grp == siliconGrp::snsr)
	{
	  // every sensor has a different label
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000  + stave*10000 + sensor*10;
	  return label_base;
	}
      if(si_grp == siliconGrp::stv)
	{
	  // layer and stave, assign all sensors to the stave number
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000 + stave*10000;
	  return label_base;
	}
      if(si_grp == siliconGrp::brrl)
	// layer only, assign all sensors to sensor 0 
	label_base += layer*1000000 + 0;
      return label_base;
    }
  else if(layer > 6 && layer < 55)
    {
      if(tpc_grp == tpcGrp::htst)
	{
	  // want every hitset (layer, sector, side) to have a separate label
	  // each group of 12 subsurfaces (sensors) is in a single hitset
	  int hitset = sensor/12; // 0-11 on side 0, 12-23 on side 1
	  label_base += layer*1000000 + hitset*10000;
	  return label_base;
	}
      if(tpc_grp == tpcGrp::sctr)
	{
	  // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
	  int side = sensor / 144; // 0-143 on side 0, 144-287 on side 1
	  int sector = (sensor - side *144) / 12; 
	  // for a given layer there are only 12 sectors x 2 sides
	  // The following gives the sectors in the inner, mid, outer regions unique group labels
	  int region = getTpcRegion(layer);  // inner, mid, outer
	  label_base += 7*1000000 + (region * 24 + side*12 + sector) *10000; 
	  return label_base;
	}
      if(tpc_grp == tpcGrp::tp)
	{
	  // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
	  label_base += 7*1000000 + 0;
	  return label_base;
	}
    }
  else
    {
      if(mms_grp == mmsGrp::tl)
	{
	  // every tile has different label
	  int tile = sensor;
	  label_base += layer*1000000 + tile*10000+sensor*10;
	  return label_base;
	}
      if(mms_grp == mmsGrp::mm)
	{
	  // assign layer 55 and tile 0 to all
	  label_base += 55*1000000 + 0;	  
	  return label_base;
	}
    }

  return -1;
}

// this method to be replaced by calls to TrackFitUtils
void HelicalFitter::getTrackletClusters(TrackSeed *tracklet, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  getTrackletClusterList(tracklet, cluskey_vec);
  // store cluster global positions in a vector
  TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   
}

void HelicalFitter::getTrackletClusterList(TrackSeed *tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = _cluster_map->findCluster(key);
      if(!cluster)
	{
	  std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}	  
      
      /// Make a safety check for clusters that couldn't be attached to a surface
      auto surf = _tGeometry->maps().getSurface(key, cluster);
      if(!surf)
	{ continue; }
      
      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}

std::vector<float> HelicalFitter::fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  return TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
}

Acts::Vector3 HelicalFitter::getClusterError(TrkrCluster *cluster, TrkrDefs::cluskey cluskey, Acts::Vector3& global)
{
  Acts::Vector3 clus_sigma(0,0,0);

  if(_cluster_version==3)
    {
      clus_sigma(2) = cluster->getZError();
      clus_sigma(0) = cluster->getRPhiError() / sqrt(2);
      clus_sigma(1) = cluster->getRPhiError() / sqrt(2);
    }
  else if(_cluster_version==4)
    {
      double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
      auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster,clusRadius,cluskey);
      float exy2 = para_errors.first * Acts::UnitConstants::cm2;
      float ez2 = para_errors.second * Acts::UnitConstants::cm2;
      clus_sigma(2) = sqrt(ez2);
      clus_sigma(0) = sqrt(exy2 / 2.0);
      clus_sigma(1) = sqrt(exy2 / 2.0);
    }
  else if(_cluster_version == 5)
    {
      double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
      TrkrClusterv5* clusterv5 = dynamic_cast<TrkrClusterv5*>(cluster);
      auto para_errors = _ClusErrPara.get_clusterv5_modified_error(clusterv5,clusRadius,cluskey);
      double phierror = sqrt(para_errors.first);
      double zerror = sqrt(para_errors.second);
      clus_sigma(2) = zerror;
      clus_sigma(0) = phierror/ sqrt(2);
      clus_sigma(1) = phierror/ sqrt(2);
    }

  return clus_sigma; 
}

void HelicalFitter::getLocalDerivativesX(Surface surf, Acts::Vector3 fitpoint, Acts::Vector3& fitpoint_local, std::vector<float>& fitpars, float lcl_derivative[5])
{
  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];
  float x = fitpoint_local(0);
  float y = fitpoint_local(2);

  // fitpoint is in global coords, fitpoint local in local coords
  // local x is unaffected by the z fit 
  // we only need to consider the circle fit paramaters
  // Do these numerically

  // dx/dradius
  // increasing R changes both local x and local y very little!
  float dx_dr = 0;  
  /*
    OR: change R and transform fitpoint to local coords
  */

  // dx/dx0
  // changing x0 changes global x by the same amount
  float dx0 = 0.01;
  Acts::Vector3 fitpoint_now(fitpoint(0)+dx0, fitpoint(1), fitpoint(2));
  Acts::Vector3 fitpoint_local_now = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (fitpoint_now *  Acts::UnitConstants::cm);
  fitpoint_local_now /= Acts::UnitConstants::cm;
  Acts::Vector3 dres = fitpoint_local_now - fitpoint_local;
  float dx_dx0 = dres(0) / dx0;  

  float dy0 = 0.01;
  // changing y0 changes global y by the same amount
  fitpoint_now(0) = fitpoint(0);
  fitpoint_now(1) = fitpoint(1)+dy0;
  fitpoint_now(2) = fitpoint(2);
  fitpoint_local_now = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (fitpoint_now *  Acts::UnitConstants::cm);
  fitpoint_local_now /= Acts::UnitConstants::cm;
  dres = fitpoint_local_now - fitpoint_local;
  float dx_dy0 = dres(0) / dy0;  

  if(Verbosity() > 1) {
    std::cout << "    x " << x << " y " << y << " x0 " << x0 << " y0 " << y0 << " R " << radius << std::endl;
    std::cout << "    LclDerivsX: dx_dx0 " << dx_dx0 << " dx_dy0 " << dx_dy0 << " dx_dr " << dx_dr << std::endl;
  }

  for(int i=0;i<5;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[0] = dx_dx0; 
  lcl_derivative[1] = dx_dy0;
  lcl_derivative[3] = dx_dr;
}

void HelicalFitter::getLocalDerivativesY(Acts::Vector3& pca, std::vector<float>& fitpars, float lcl_derivative[5])
{
  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];
  float x = pca(0);
  float y = pca(1);
  float dr = 0.2;
  float dx0 = 0.05;

  // dy/dradius
  float phi = atan2(y-y0, x-x0);
  float dy = dr * sin(phi);
  float dy_dr = dy/dr;

  // dy/dy0
  float dy_dy0 = 1.0;

  // dy/dx0
  // y = y0 + sqrt(pow(radius, 2) + pow(x-x0, 2))  
  float dyx0 = sqrt(pow(radius, 2) + pow(x-x0+dx0, 2)) 
                      - sqrt(pow(radius, 2) + pow(x-x0, 2));
  float dy_dx0 = dyx0/dx0;

  if(Verbosity() > 1) { 
    std::cout << "    x " << x << " y " << y << " x0 " << x0 << " y0 " << y0 << " R " << radius << std::endl;
    std::cout << "    LclDerivsY: dy_dy0 " << dy_dy0 << " dy_dx0 " << dy_dx0 << " dy_dr " << dy_dr << std::endl; 
  }
  
  for(int i=0;i<NLC;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[0] = dy_dx0;
  lcl_derivative[1] = dy_dy0;
  lcl_derivative[3] = dy_dr;
}

void HelicalFitter::getLocalDerivativesZ(Acts::Vector3& fitpoint,  std::vector<float>& fitpars, float lcl_derivative[5])
{
  // the local coord corresponding to z is local-y.
  // changes in z global translate directly to y-local

  float zslope = fitpars[3];
  float cluster_radius = sqrt(fitpoint(0)*fitpoint(0)+fitpoint(1)*fitpoint(1));
  float z0 = fitpars[4];

  // z = z0 + zslope * cluster_radius 
  float dz_dz0 = 1.0;
  float dz_dzslope = cluster_radius;
  //  float dz_dradius = zslope;
  if(Verbosity() > 1) {
    std::cout << "    x " << fitpoint(0)<<" y "<<fitpoint(1)<<" z "<<fitpoint(2)<< " z0 " <<z0 << " zslope " <<zslope <<std::endl;  
    std::cout << "    LclDerivsZ: dz_dzslope " << dz_dzslope << " dz_dz0 " << dz_dz0 
      //<< " dz_dradius " << dz_dradius 
	      << std::endl;
  }
  
  for(int i=0;i<NLC;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[2] = dz_dz0;   // dz/dz_0
  //lcl_derivative[3] = dz_dradius;   // dz/dradius
  lcl_derivative[4] = dz_dzslope;   // dz/dz_slope
}

 void HelicalFitter::getGlobalDerivativesX( std::vector<Acts::Vector3> angleDerivs, std::vector<Acts::Vector3> translDerivs, float glbl_derivative[], unsigned int layer)
{
  // local-x:  relevant global pars are alpha, beta, gamma, dx, dy (ipar 0,1,2,3)
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer))
    {
      //glbl_derivative[3] = 1.0;  // optimize dx
      glbl_derivative[0] = angleDerivs[0](0);  // dx/dalpha
      glbl_derivative[1] = angleDerivs[1](0);  // dx/dbeta
      glbl_derivative[2] = angleDerivs[2](0);  // dx/dgamma

      glbl_derivative[3] = translDerivs[0](0);  // dx/dx
      glbl_derivative[4] = translDerivs[1](0);  // dx/dy = 0
      glbl_derivative[5] = translDerivs[2](0);  // dx/dz = 0

      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    }
}

void HelicalFitter::getGlobalDerivativesY( std::vector<Acts::Vector3> angleDerivs,  std::vector<Acts::Vector3> translDerivs, float glbl_derivative[], unsigned int layer)
{
  // y - relevant global pars are alpha, beta, gamma, dy (ipar 0,1,2,4)
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer)) 
    {
      //glbl_derivative[4] = 1.0; // optimize dy
      glbl_derivative[0] = angleDerivs[0](1);   // dy/dalpha
      glbl_derivative[1] = angleDerivs[1](1);   // dy/dbeta
      glbl_derivative[2] = angleDerivs[2](1);   // dy/dgamma

      glbl_derivative[3] = translDerivs[0](1);  // dy/dx
      glbl_derivative[4] = translDerivs[1](1);  // dy/dy
      glbl_derivative[5] = translDerivs[2](1);  // dy/dz
      
      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    }
}

void HelicalFitter::getGlobalDerivativesZ( std::vector<Acts::Vector3> angleDerivs,  std::vector<Acts::Vector3> translDerivs, float glbl_derivative[], unsigned int layer)
{
  // z - relevant global pars are alpha, beta, dz (ipar 0,1,5)
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer))
    { 
      //glbl_derivative[5] = 1.0;  // optimize dz
      glbl_derivative[0] = angleDerivs[0](2);  // dz/dalpha
      glbl_derivative[1] = angleDerivs[1](2);  // dz/dbeta
      glbl_derivative[2] = angleDerivs[2](2);  // dz/dgamma

      glbl_derivative[3] = translDerivs[0](2);  // dz/dx = 0
      glbl_derivative[4] = translDerivs[1](2);  // dz/dy = 0
      glbl_derivative[5] = translDerivs[2](2);  // dz/dz

      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    } 
}

void HelicalFitter::printBuffers(int index, Acts::Vector2 residual, Acts::Vector3 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
{
  std::cout << " float buffer: " << " residual " << "  " << residual(index);
  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_deriv["<< il << "] " << lcl_derivative[il] << "  ";  }
  std::cout  << " sigma " << "  " << clus_sigma(index) << "  ";
  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  ";  }
  std::cout << " int buffer: " << " 0 " << "  ";
  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_label["<< il << "] " << il << "  ";  }
  std::cout << " 0 " << "  ";
  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  ";  }
  std::cout << " end of meas " << std::endl;		    
}

unsigned int HelicalFitter::addSiliconClusters(std::vector<float>& fitpars, std::vector<Acts::Vector3>& global_vec,  std::vector<TrkrDefs::cluskey>& cluskey_vec)
{

  return TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
}

 bool HelicalFitter::is_layer_fixed(unsigned int layer)
 {
  bool ret = false;
  auto it = fixed_layers.find(layer);
  if(it != fixed_layers.end()) 
    ret = true;

  return ret;
 }

 void HelicalFitter::set_layer_fixed(unsigned int layer)
 {
   fixed_layers.insert(layer);
 }
 
bool HelicalFitter::is_layer_param_fixed(unsigned int layer, unsigned int param)
 {
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = fixed_layer_params.find(pair);
  if(it != fixed_layer_params.end()) 
    ret = true;

  return ret;
 }

void HelicalFitter::set_layer_param_fixed(unsigned int layer, unsigned int param)
 {
   std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
   fixed_layer_params.insert(pair);
 }

void HelicalFitter::correctTpcGlobalPositions(std::vector<Acts::Vector3> global_vec,  std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  for(unsigned int iclus=0;iclus<cluskey_vec.size();++iclus)
    {
      auto cluskey = cluskey_vec[iclus];
      auto global = global_vec[iclus];
      const unsigned int trkrId = TrkrDefs::getTrkrId(cluskey);	  
      if(trkrId == TrkrDefs::tpcId) 
	{  
	  // have to add corrections for TPC clusters after transformation to global
	  int crossing = 0;  // for now
	  makeTpcGlobalCorrections(cluskey, crossing, global); 
	  global_vec[iclus] = global;
	}
    }
}


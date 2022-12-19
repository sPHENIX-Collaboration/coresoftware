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
	 << " TPC track map size " << _track_map_tpc->size() 
	 << " Silicon track map size "  << _track_map_silicon->size() 
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
      getTrackletClusters(tracklet, global_vec, cluskey_vec);   // store cluster corrected global positions in a vector
      std::vector<float> fitpars =  fitClusters(global_vec, cluskey_vec);       // do helical fit
      if(fitpars.size() == 0) continue;  // discard this track, not enough clusters to fit

      if(Verbosity() > 0)  
	{ std::cout << " Track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
		 << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }

      // if a full track is requested, get the silicon clusters too and refit
      if(fittpc && fitfulltrack)
	{
	  // this associates silicon clusters and adds them to the vectors
	  unsigned int nsilicon = addSiliconClusters(fitpars, global_vec, cluskey_vec);
	  if(nsilicon < 3) continue;  // discard this TPC seed, did not get a good match to silicon

	  // fit the full track now
	  fitpars.clear();
	  fitpars =  fitClusters(global_vec, cluskey_vec);       // do helical fit
	  if(fitpars.size() == 0) continue;  // discard this track, fit failed

	  if(Verbosity() > 0)  
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
	  
	  // PCA of helix to cluster global position
	  Acts::Vector3 pca = get_helix_pca(fitpars, global);
	  if(Verbosity() > 0) {std::cout << "    cluster position " << global(0) << " " << global(1) << " " << global(2) 
					 << " pca " << pca(0) << " " << pca(1) << " " << pca(2)  << std::endl;}
	    
	  // capture residuals in the form of (data - fit)
	  auto residual = global - pca;
	  // need standard deviation of measurements
	  Acts::Vector3 clus_sigma = getClusterError(cluster, cluskey, global);
	  if(isnan(clus_sigma(0)) || isnan(clus_sigma(1)) || isnan(clus_sigma(2)))  { continue; }

	  Surface surf = _tGeometry->maps().getSurface(cluskey, cluster);
	  int glbl_label[NGL];
	  getGlobalLabels(surf, glbl_label);  // these depend on the sensor grouping

	  float lcl_derivative[NLC];
	  float glbl_derivative[NGL];
	  unsigned int crossing = 0;
	  // The angleDerivs dimensions are [alpha/beta/gamma](x/y/z)
	  std::vector<Acts::Vector3> angleDerivs = getDerivativesAlignmentAngles(global, cluskey, cluster, surf, crossing); 

	  // Add the measurement separately for each coordinate direction to Mille
	  // set the derivatives non-zero only for parameters we want to be optimized
	  unsigned int layer = TrkrDefs::getLayer(cluskey_vec[ivec]);	  
	  getLocalDerivativesX(pca, fitpars, lcl_derivative);
	  getGlobalDerivativesX(angleDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3)
	    { std::cout << "layer " << layer << " X buffers:" << std::endl; printBuffers(0, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if( !isnan(residual(0)) && clus_sigma(0) < 1.0)  // discards crazy clusters
	    { _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(0), clus_sigma(0));}

	  getLocalDerivativesY(pca, fitpars, lcl_derivative);
	  getGlobalDerivativesY(angleDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3) 
	    { std::cout  << "layer " << layer << " Y buffers:" << std::endl; printBuffers(1, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if( !isnan(residual(1)) && clus_sigma(1) < 1.0)  // discards crazy clusters
	    {_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(1), clus_sigma(1));}

	  getLocalDerivativesZ(pca, lcl_derivative);
	  getGlobalDerivativesZ(angleDerivs, glbl_derivative, layer);
	  if(Verbosity() > 3) 
	    { std::cout  << "layer " << layer << " Z buffers:" << std::endl; printBuffers(2, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); }
	  if(!isnan(residual(2)) && clus_sigma(2) < 1.0)
	    {_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(2), clus_sigma(2));}
	}

      // close out this track
      _mille->end();
      
    }  // end loop over tracks
  
  return Fun4AllReturnCodes::EVENT_OK;
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

std::vector<Acts::Vector3> HelicalFitter::getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing)
{
  // The value of global is from the geocontext transformation
  // we add to that transformation a small rotation around the relevant axis in the surface frame

  std::vector<Acts::Vector3> derivs_vector;

  // get the transformation from the geocontext
  Acts::Transform3 transform = surface->transform(_tGeometry->geometry().getGeoContext());

  // Make an additional transform that applies a small rotation angle in the surface frame
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
	  Acts::Transform3 overallTransformation = transform * perturbationTransformation;
	  
	  // transform the cluster local position to global coords with this additional rotation added
	  auto x = cluster->getLocalX() * 10;  // mm 
	  auto y = cluster->getLocalY() * 10;	  
	  if(trkrId == TrkrDefs::tpcId) { y = convertTimeToZ(cluster_key, cluster); }
	  
	  Eigen::Vector3d clusterLocalPosition (x,y,0);  // follows the convention for the acts transform of local = (x,z,y)
	  Eigen::Vector3d finalCoords = overallTransformation*clusterLocalPosition;  // result in mm
	  finalCoords /= 10.0; // convert mm back to cm
	  
	  // have to add corrections for TPC clusters after transformation to global
	  if(trkrId == TrkrDefs::tpcId) {  makeTpcGlobalCorrections(cluster_key, crossing, global); }

	  // note that global cancels out here	  	  
	  if(ip == 0)
	    {
	      keeper(0) = (finalCoords(0) - global(0));
	      keeper(1) = (finalCoords(1) - global(1));
	      keeper(2) = (finalCoords(2) - global(2));
	    }
	  else
	    {
	      keeper(0) -= (finalCoords(0) - global(0));
	      keeper(1) -= (finalCoords(1) - global(1));
	      keeper(2) -= (finalCoords(2) - global(2));
	    }

	  if(Verbosity() > 1)
	    {
	      std::cout << "        finalCoords(0) " << finalCoords(0) << " global(0) " << global(0) << " finalCoords(1) " << finalCoords(1) 
			<< " global(1) " << global(1) << " finalCoords(2) " << finalCoords(2) << " global(2) " << global(2) << std::endl;
	      std::cout  << "        keeper now:  keeper(0) " << keeper(0) << " keeper(1) " << keeper(1) << " keeper(2) " << keeper(2) << std::endl;
	    }
	}

      // derivs vector contains:
      //   (dx/dalpha,     dy/dalpha,     dz/dalpha)     (for iangle = 0)
      //   (dx/dbeta,        dy/dbeta,        dz/dbeta)        (for iangle = 1) 
      //   (dx/dgamma, dy/dgamma, dz/dgamma) (for iangle = 2)
 
     // Average the changes to get the estimate of the derivative      
      derivs(0) = keeper(0) / (2.0 * fabs(theseAngles[iangle]));
      if( isnan(derivs(0)) ) { derivs(0) = 0; }
      derivs(1) = keeper(1) / (2.0 * fabs(theseAngles[iangle]));
      if( isnan(derivs(1)) ) { derivs(1) = 0; }
      derivs(2) = keeper(2) / (2.0 * fabs(theseAngles[iangle]));
      if( isnan(derivs(2)) ) { derivs(2) = 0; }
      derivs_vector.push_back(derivs);

      if(Verbosity() > 1) { std::cout << "        derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl; }
    }
  
  return derivs_vector;
}

  Acts::Transform3 HelicalFitter::makePerturbationTransformation(Acts::Vector3 angles)
  {
    // Note: Here beta is apllied to the z axis and gamma is applied to the y axis because the geocontext transform 
    // will flip those axes when transforming to global coordinates
    Eigen::AngleAxisd alpha(angles(0), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd beta(angles(2), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd gamma(angles(1), Eigen::Vector3d::UnitZ());
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

void HelicalFitter::getTrackletClusters(TrackSeed *tracklet, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec)
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
      
      Acts::Vector3 global  = _tGeometry->getGlobalPosition(key, cluster);	  
      
      const unsigned int trkrId = TrkrDefs::getTrkrId(key);	  
      
      // have to add corrections for TPC clusters after transformation to global
      if(trkrId == TrkrDefs::tpcId) 
	{  
	  int crossing = 0;  // for now
	  makeTpcGlobalCorrections(key, crossing, global); 
	}
      
      // add the global positions to a vector to give to the helical fitter
      global_vec.push_back(global);
      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}

std::vector<float> HelicalFitter::fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec)
{
     std::vector<float> fitpars;

      // make the helical fit using TrackFitUtils
      if(global_vec.size() < 3)  
	if(Verbosity() > 0) {  std::cout << " track has too few clusters for circle fit, skip it" << std::endl; return fitpars; }
      std::tuple<double, double, double> circle_fit_pars = TrackFitUtils::circle_fit_by_taubin(global_vec);

      // It is problematic that the large errors on the INTT strip z values are not allowed for - drop the INTT from the z line fit
      std::vector<Acts::Vector3> global_vec_noINTT;
      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey_vec[ivec]);
	  if(trkrid != TrkrDefs::inttId) { global_vec_noINTT.push_back(global_vec[ivec]); }
	}      
      if(global_vec_noINTT.size() < 3) 
	if(Verbosity() > 0) { std::cout << " track has too few non-INTT clusters for z fit, skip it" << std::endl; return fitpars; }
     std::tuple<double,double> line_fit_pars = TrackFitUtils::line_fit(global_vec_noINTT);

     fitpars.push_back( std::get<0>(circle_fit_pars));
     fitpars.push_back( std::get<1>(circle_fit_pars));
     fitpars.push_back( std::get<2>(circle_fit_pars));
     fitpars.push_back( std::get<0>(line_fit_pars));
     fitpars.push_back( std::get<1>(line_fit_pars));

     return fitpars; 
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
  return clus_sigma; 
}

void HelicalFitter::getLocalDerivativesX(Acts::Vector3& pca, std::vector<float>& fitpars, float lcl_derivative[5])
{
  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];

  float x_x0 = 1.0;
  // dx/dradius = radius / (R^2 - (y-y0)^2)^1/2  , sign of (x-x0)/fabs(x-x0)
  float x_radius = ( pca(0)-x0 ) / fabs( pca(0)-x0 ) * radius / sqrt(radius*radius - (pca(1) - y0 ) *  (pca(1) - y0 ) );
  if(Verbosity() > 0) {std::cout << "    x derivatives: x_x0 " << x_x0 << " x_radius " << x_radius << std::endl;}
  if(isnan(x_radius)) {x_radius = 0.0; }
  // dx/dy0 = (y-y0) / (R^2 - (y-y0)^2)^1/2 
  float x_y0 = (pca(1) - y0) / sqrt(radius*radius - (pca(1) - y0) * (pca(1) - y0));
  if(isnan(x_y0)) {x_y0 = 0.0; }
  
  for(int i=0;i<5;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[0] = x_x0; 
  lcl_derivative[1] = x_y0;
  lcl_derivative[3] = x_radius;
}

void HelicalFitter::getLocalDerivativesY(Acts::Vector3& pca, std::vector<float>& fitpars, float lcl_derivative[5])
{
  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];

  float y_y0 = 1.0;
  // dy/dradius = sign * radius / (R^2 - (x-x0)^2)^1/2 , sign of (y-y0)/fabs(y-y0)
  float y_radius =  ( pca(1)-y0 ) / fabs( pca(1)-y0 ) * radius / sqrt(radius*radius - (pca(0) - x0 ) * (pca(0) - x0 ) );
  if(Verbosity() > 0) { std::cout << "    y derivatives: y_y0 " << y_y0 << " y_radius " << y_radius << std::endl; }
  if(isnan(y_radius)) {y_radius = 0.0; }
  // dy/dx0 = (x-x0) / (R^2 - (x-x0)^2)^1/2 
  float y_x0 = (pca(0) - x0) / sqrt(radius*radius - (pca(0) - x0) * (pca(0) - x0));
  if(isnan(y_x0)) {y_x0 = 0.0; }
  
  for(int i=0;i<NLC;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[0] = y_x0;
  lcl_derivative[1] = y_y0;
  lcl_derivative[3] = y_radius;
}

void HelicalFitter::getLocalDerivativesZ(Acts::Vector3& global, float lcl_derivative[5])
{
  // z = z0 + zslope * cluster_radius 
  float cluster_radius = sqrt(global(0)*global(0)+global(1)*global(1));
  float z_zslope = cluster_radius;
  float z_z0 = 1.0;
  if(Verbosity() > 0) {std::cout << "    z derivatives: z_zslope " << z_zslope << " z_z0 " << z_z0 << std::endl;}
  
  for(int i=0;i<NLC;++i) {lcl_derivative[i] = 0.0;}
  lcl_derivative[2] = z_z0;
  lcl_derivative[4] = z_zslope;
}

void HelicalFitter::getGlobalDerivativesX( std::vector<Acts::Vector3> angleDerivs, float glbl_derivative[], unsigned int layer)
{
  // x - relevant global pars are alpha, beta, gamma, dx (ipar 0,1,2,3), relevant local pars are x_x0, x_radius
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer))
    {
      glbl_derivative[3] = 1.0;  // optimize dx
      glbl_derivative[0] = angleDerivs[0](0);  // dx/dalpha
      glbl_derivative[1] = angleDerivs[1](0);  // dx/dbeta
      glbl_derivative[2] = angleDerivs[2](0);  // dx/dgamma

      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    }
}

void HelicalFitter::getGlobalDerivativesY( std::vector<Acts::Vector3> angleDerivs, float glbl_derivative[], unsigned int layer)
{
  // y - relevant global pars are alpha, beta, gamma, dy (ipar 0,1,2,4)
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer)) 
    {
      glbl_derivative[4] = 1.0; // optimize dy
      glbl_derivative[0] = angleDerivs[0](1);   // dy/dalpha
      glbl_derivative[1] = angleDerivs[1](1);   // dy/dbeta
      glbl_derivative[2] = angleDerivs[2](1);   // dy/dgamma
      
      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    }
}

void HelicalFitter::getGlobalDerivativesZ( std::vector<Acts::Vector3> angleDerivs, float glbl_derivative[], unsigned int layer)
{
  // z - relevant global pars are alpha, beta, dz (ipar 0,1,5)
  for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
  if(!is_layer_fixed(layer))
    { 
      glbl_derivative[5] = 1.0;  // optimize dz
      glbl_derivative[0] = angleDerivs[0](2);  // dz/dalpha
      glbl_derivative[1] = angleDerivs[1](2);  // dz/dbeta
      glbl_derivative[2] = angleDerivs[2](2);  // dz/dgamma

      for(int i=0; i< NGL; ++i)
	{
	  if(is_layer_param_fixed(layer,i))
	    {glbl_derivative[i] = 0.0;}
	}
    } 
}

void HelicalFitter::printBuffers(int index, Acts::Vector3 residual, Acts::Vector3 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
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
  // project the fit of the TPC clusters to each silicon layer, and find the nearest silicon cluster
  // iterate over the cluster map and find silicon clusters that match this track fit

  unsigned int nsilicon = 0;

  // We want the best match in each layer
  std::vector<float> best_layer_dca;
  best_layer_dca.assign(7, 999.0);
  std::vector<TrkrDefs::cluskey> best_layer_cluskey;
  best_layer_cluskey.assign(7, 0);

  for(const auto& hitsetkey:_cluster_map->getHitSetKeys())
    {
      auto range = _cluster_map->getClusters(hitsetkey);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  TrkrDefs::cluskey cluskey = clusIter->first;
	  unsigned int layer = TrkrDefs::getLayer(cluskey);
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
	  
	  if(trkrid != TrkrDefs::mvtxId && trkrid != TrkrDefs::inttId)  continue;
	  
	  TrkrCluster* cluster = clusIter->second;
	  auto global = _tGeometry->getGlobalPosition(cluskey, cluster);

	  Acts::Vector3 pca = get_helix_pca(fitpars, global);
	  float dca = (pca - global).norm();
	  if(trkrid == TrkrDefs::inttId || trkrid == TrkrDefs::mvtxId)
	    {
	      Acts::Vector2 global_xy(global(0), global(1));
	      Acts::Vector2 pca_xy(pca(0), pca(1));
	      Acts::Vector2 pca_xy_residual = pca_xy - global_xy;
	      dca = pca_xy_residual.norm();
	    }

	  if(dca < best_layer_dca[layer])
	    {
	      best_layer_dca[layer] = dca;
	      best_layer_cluskey[layer] = cluskey;
	    }
	}  // end cluster iteration
    } // end hitsetkey iteration

  for(unsigned int layer = 0; layer < 7; ++layer)
    {
      if(best_layer_dca[layer] < dca_cut)
	{
	  cluskey_vec.push_back(best_layer_cluskey[layer]);
	  auto clus =  _cluster_map->findCluster(best_layer_cluskey[layer]);
	  auto global = _tGeometry->getGlobalPosition(best_layer_cluskey[layer], clus);
	  global_vec.push_back(global);
	  nsilicon++;
	  if(Verbosity() > 0) std::cout << "   add cluster in layer " << layer << " with cluskey " << best_layer_cluskey[layer] << " and dca " << best_layer_dca[layer] 
		    << " nsilicon " << nsilicon << std::endl;
	}
    }

  return nsilicon;
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

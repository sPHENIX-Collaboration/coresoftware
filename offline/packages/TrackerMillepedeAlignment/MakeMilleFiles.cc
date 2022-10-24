#include "MakeMilleFiles.h"

#include "Mille.h"

/// Tracking includes

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrClusterv3.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/ActsTransformations.h>

#include <trackbase/TpcDefs.h>               // for side

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

#include <TF1.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <fstream>

//____________________________________________________________________________..
MakeMilleFiles::MakeMilleFiles(const std::string &name)
  : SubsysReco(name)
  , _mille(nullptr)
{}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  //  _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  _mille = new Mille(data_outfilename.c_str()); 

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename.c_str());
  steering_file << data_outfilename.c_str() << std::endl;
  steering_file.close();

  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are silicon " << si_group << " tpc " << tpc_group << " mms " << mms_group << std::endl; 


  return ret;
}

//____________________________________________________________________________..
int MakeMilleFiles::process_event(PHCompositeNode */*topNode*/)
{
  // Outline:
  //
  // loop over tracks 
  //   Make any track cuts here to skip undesirable tracks (maybe low pT?)
  //   loop over measurements for each track
  //      for each measurement
  //         Get measurement value and error (global, what to use for error?)
  //         Get track state at measurement and determine residual, error is measurement error
  //         Determine # of local parameters
  //         Determine # of global parameters and their labels
  //         Loop over local parameters
  //            Determine derivative of residual wrt to this parameter
  //            Add to array of local derivatives, array index is used as local parameter label
  //         Loop over global parameters
  //            Determine derivatives of residual wrt this parameter
  //            Add to array of global derivatives, add global par label to array of labels
  //          Call _mille->mille() with arguments:
  //              #local pars
  //              array of local derivatives
  //              #global pars
  //              array of global derivatives
  //              array of integer global par labels
  //              residual value (float) z = measurement - track state
  //              sigma of measurement
  //   After processing all measurements for this track, call _mille->end() to add buffer to file and reset buffer
  // After all tracks are processed, file is closed when Mille destructor is called


  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;

  // loop over the tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      auto track = phtrk_iter->second;
      auto crossing = track->get_silicon_seed()->get_crossing(); 
      
      if (Verbosity() > 0)
	{
	  std::cout << std::endl << __LINE__   << ": Processing track itrack: " << phtrk_iter->first << ": nhits: " << track-> size_cluster_keys()
		    << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
	}

      // Make any desired track cuts here
      // Maybe set a lower pT limit - low pT tracks are not very sensitive to alignment
     
      std::map<TrkrDefs::cluskey, TrkrCluster*> all_clusters;
      for( const auto& seed: { track->get_silicon_seed(), track->get_tpc_seed() } )
	{
      // Get all clusters for this track from the seeds
	  for (auto key_iter = seed->begin_cluster_keys();
	       key_iter != seed->end_cluster_keys();
	       ++key_iter)
	    {
	      TrkrDefs::cluskey cluster_key = *key_iter;
	      auto cluster = _cluster_map->findCluster(cluster_key);
	      if(!cluster) { continue;}
	      all_clusters.insert(std::make_pair(cluster_key, cluster));
	    }

	  for(auto clus_iter : all_clusters)
	    {
	      TrkrDefs::cluskey cluster_key = clus_iter.first;
	      auto cluster = clus_iter.second;
 	      
	      // we want the global cluster positions
	      Acts::Vector3 global  = _tGeometry->getGlobalPosition(cluster_key, cluster);
	      
	      // TPC clusters need distortion corrections, silicon and MM's clusters do not
	      unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
	      if(trkrId == TrkrDefs::tpcId) { makeTpcGlobalCorrections(cluster_key, crossing, global); }
	      
	      // we have our global cluster position, corrected if necessary. Each component is three measurements, x, y and z
	      // we need to find the residual and its derivative wrt parameters in each coordinate direction
	      
	      // find track state that is the closest to cluster 
	      auto state_iter =  getStateIter(global, track);

	      // Use the PCA of the track to the measurement to get the residual
	      Acts::Vector3 pca = getPCALinePoint(global, state_iter->second);
	      Acts::Vector3 residual = global - pca;
	      
	      if(Verbosity() > 1)
		{
		  std::cout << " cluster global " << global(0) << "  " << global(1) << "  " << global(2) << " track PCA " << pca(0) << "  " << pca(1) << "  " << pca(2) 
			    << " residual " << residual(0) << "  " << residual(1) << "  " << residual(2)  << std::endl;
		}
	      
	      // need standard deviation of measurements
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
		  auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster,clusRadius,cluster_key);
		  float exy2 = para_errors.first * Acts::UnitConstants::cm2;
		  float ez2 = para_errors.second * Acts::UnitConstants::cm2;
		  clus_sigma(2) = sqrt(ez2);
		  clus_sigma(0) = sqrt(exy2 / 2.0);
		  clus_sigma(1) = sqrt(exy2 / 2.0);
		}
	      
	      // Get the surface for this cluster
	      Surface surf = _tGeometry->maps().getSurface(cluster_key, cluster);

	      // if this is a TPC cluster, check that the corrections did not change the surface
	      // ***** Do we really want to do this? We are aligning the readout pads
	      // ***** BUT: the tracker uses the transformation for the new surface, not the readout surface
	      if(trkrId == TrkrDefs::tpcId)
		{
		  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
		  TrkrDefs::subsurfkey new_subsurfkey = 0;    
		  surf = _tGeometry->get_tpc_surface_from_coords(hitsetkey,  global, new_subsurfkey);
		}
	      if(!surf)  { continue; }

	      // For now, ignore local pars
	      static const int NLC = 5;
	      float lcl_derivative[NLC] = {0,0,0,0,0};
	      	      
	      // The global alignment parameters are given initial values of zero by default, we do not specify them
	      // We identify the global alignment parameters for this surface
	      Acts::GeometryIdentifier id = surf->geometryId();
	      int label_base = getLabelBase(id);   // This value depends on how the surfaces are grouped
	      
	      static const int NGL = 6;
	      int glbl_label[NGL];
	      for(int i=0;i<NGL;++i) 
		{
		  glbl_label[i] = label_base + i;
		  if(Verbosity() > 1) { std::cout << "  glbl " << i << " label " << glbl_label[i] << " "; }
		}
	      if(Verbosity() > 1) { std::cout << std::endl; }

	      // Add the measurement separately for each coordinate direction to Mille
	      // set the derivatives non-zero only for parameters we want to be optimized
	      float glbl_derivative[NGL];
	      // The angleDerivs dimensions are [alpha/beta/gamma](x/y/z)
	      std::vector<Acts::Vector3> angleDerivs = getDerivativesAlignmentAngles(global, cluster_key, cluster, surf, crossing); 

	      // x - relevant global pars are alpha, beta, gamma, dx (ipar 0,1,2,3)
	      for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
	      glbl_derivative[3] = 1.0;  // optimize dx
	      glbl_derivative[0] = angleDerivs[0](0);  // dx/dalpha
	      glbl_derivative[1] = angleDerivs[1](0);  // dx/dbeta
	      glbl_derivative[2] = angleDerivs[2](0);  // dx/dgamma
	      if(clus_sigma(0) < 1.0)  // discards crazy clusters
		{ _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(0), clus_sigma(0));}

	      if(Verbosity() > 3)
		{
		  std::cout << " X:  float buffer: " << " residual " << "  " << residual(0);
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_deriv["<< il << "] " << lcl_derivative[il] << "  ";  }
		  std::cout  << " sigma " << "  " << clus_sigma(0) << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  ";  }
		  std::cout << " X:  int buffer: " << " 0 " << "  ";
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_label["<< il << "] " << il << "  ";  }
		  std::cout << " 0 " << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  ";  }
		  std::cout << " end of X meas " << std::endl;		    
		}
	      
	      // y - relevant global pars are alpha, beta, gamma, dy (ipar 0,1,2,4)
	      for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
	      glbl_derivative[4] = 1.0; // optimize dy
	      glbl_derivative[0] = angleDerivs[0](1);   // dy/dalpha
	      glbl_derivative[1] = angleDerivs[1](1);   // dy/dbeta
	      glbl_derivative[2] = angleDerivs[2](1);   // dy/dgamma
	      if(clus_sigma(1) < 1.0)  // discards crazy clusters
		{_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(1), clus_sigma(1));}
	      
	      if(Verbosity() > 3) 
		{ 
		  std::cout << " Y:  float buffer: " << " residual " << "  " << residual(1);
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_deriv["<< il << "] " << lcl_derivative[il] << "  ";  }
		  std::cout  << " sigma " << "  " << clus_sigma(1) << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  ";  }
		  std::cout << " Y:  int buffer: " << " 0 " << "  ";
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_label["<< il << "] " << il << "  "; }
		  std::cout << " 0 " << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  "; }
		  std::cout << " end of Y meas " << std::endl;		    
		}
	      
	      // z - relevant global pars are alpha, beta, dz (ipar 0,1,5)
	      glbl_derivative[5] = 1.0;  // optimize dz
	      glbl_derivative[0] = angleDerivs[0](2);  // dz/dalpha
	      glbl_derivative[1] = angleDerivs[1](2);  // dz/dbeta
	      glbl_derivative[2] = angleDerivs[2](2);  // dz/dgamma
	      if(clus_sigma(2) < 1.0)
		{_mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(2), clus_sigma(2));}
	      
	      if(Verbosity() > 3)
		{
		  std::cout << " Z:  float buffer: " << " residual " << "  " << residual(2);
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_deriv["<< il << "] " << lcl_derivative[il] << "  "; }
		  std::cout  << " sigma " << "  " << clus_sigma(2) << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  "; }
		  std::cout << " Z:  int buffer: " << " 0 " << "  ";
		  for (int il=0;il<NLC;++il) { if(lcl_derivative[il] != 0) std::cout << " llc_label["<< il << "] " << il << "  ";  }
		  std::cout << " 0 " << "  ";
		  for (int ig=0;ig<NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  ";  }
		  std::cout << " end of Z meas " << std::endl;		    
		}
	    }
	}
      
      // close out this track
      _mille->end();
      
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  MakeMilleFiles::GetNodes(PHCompositeNode* topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

// tpc distortion corrections
  _dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( _dcc_static )
    { 
      std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl; 
    }
  _dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  if( _dcc_average )
    { 
      std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl; 
    }
  _dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");
  if( _dcc_fluctuation )
    { 
      std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl; 
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 MakeMilleFiles::getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state)
{
  // Approximate track with a straight line consisting of the state position and the vector (px,py,pz)   

  Acts::Vector3 track_dir(state->get_px(), state->get_py(), state->get_pz());
  track_dir = track_dir / track_dir.norm(); 
  Acts::Vector3 track_base(state->get_x(), state->get_y(), state->get_z());

  // The position of the closest point on the line is:
  // track_base + projection of difference between the point and track_base on the line vector
  Acts::Vector3 pca = track_base + ( (global - track_base).dot(track_dir) ) * track_dir;

  return pca;
}

std::vector<Acts::Vector3> MakeMilleFiles::getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing)
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
	  Eigen::Vector3d finalCoords = overallTransformation*clusterLocalPosition / 10.0; // convert mm back to cm
	  
	  // have to add corrections for TPC clusters after transformation to global
	  if(trkrId == TrkrDefs::tpcId) {  makeTpcGlobalCorrections(cluster_key, crossing, global); }
	  	  
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
      // Check what the sign of this should be !!!!
      derivs(0) = keeper(0) / (2.0 * fabs(theseAngles[iangle]));
      derivs(1) = keeper(1) / (2.0 * fabs(theseAngles[iangle]));
      derivs(2) = keeper(2) / (2.0 * fabs(theseAngles[iangle]));
      derivs_vector.push_back(derivs);

      if(Verbosity() > 1)
	{
	  std::cout << "        derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl;
	}
    }
  
  return derivs_vector;
}

  Acts::Transform3 MakeMilleFiles::makePerturbationTransformation(Acts::Vector3 angles)
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

float MakeMilleFiles::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
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

void MakeMilleFiles::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
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

SvtxTrack::StateIter MakeMilleFiles::getStateIter(Acts::Vector3& global, SvtxTrack* track)
{
  float clus_radius = sqrt(global[0]*global[0]+global[1]*global[1]);
  auto state_iter = track->begin_states();
  float dr_min = -1;
  //for( auto iter = state_iter; iter != track->end_states(); ++iter )
  for( auto iter = track->begin_states(); iter != track->end_states(); ++iter )
    {
      const auto dr = std::abs( clus_radius - sqrt( iter->second->get_x()* iter->second->get_x() + iter->second->get_y()* iter->second->get_y())  );
      if( dr_min < 0 || dr < dr_min )
	{
	  state_iter = iter;
	  dr_min = dr;
	}
      else
	{ break; }
    }
  
  if(Verbosity() > 2)
    {
      std::cout << "track state for pathlength " << clus_radius << " dr_min " << dr_min << " position " 
		<< state_iter->second->get_x() << "  "
		<< state_iter->second->get_y() << "  "			
		<< state_iter->second->get_z() << std::endl;
    }  
  return state_iter;
}

int MakeMilleFiles::getLabelBase(Acts::GeometryIdentifier id)
{
  unsigned int volume = id.volume(); 
  unsigned int acts_layer = id.layer();
  unsigned int layer = base_layer_map.find(volume)->second + acts_layer / 2 -1;
  unsigned int sensor = id.sensitive() - 1;  // Acts starts at 1

  int label_base = 1;  // Mille wants to start at 1

  // decide what level of grouping we want
  if(layer < 7)
    {
      if(si_group == siliconGroup::sensor)
	{
	  // every sensor has a different label
	  label_base += layer*1000 + sensor*10;
	  return label_base;
	}
      if(si_group == siliconGroup::stave)
	{
	  // layer and stave, assign all sensors to the stave number
	  int stave = sensor / nstaves[layer];
	  label_base += layer*1000 + stave*10;
	  return label_base;
	}
      if(si_group == siliconGroup::barrel)
	// layer only, assign all sensors to sensor 0 
	label_base += layer*1000 + 0;
      return label_base;
    }
  else if(layer > 6 && layer < 55)
    {
      if(tpc_group == tpcGroup::subsurf)
	{
	  // every surface has separate label
	  label_base += layer*1000 + sensor*10;
	  return label_base;
	}
      if(tpc_group == tpcGroup::sector)
	{
	  // all tpc layers, assign layer 7 and side and sector number to all layers and subsurfaces
	  int side = sensor / 2;   // check!!!!
	  int sector = (sensor - side *144) / 12; 
	  label_base += 7*1000 + side * 1000 + sector*10; 
	  return label_base;
	}
      if(tpc_group == tpcGroup::tpc)
	{
	  // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
	  label_base += 7*1000 + 0;
	  return label_base;
	}
    }
  else
    {
      if(mms_group == mmsGroup::tile)
	{
	  // every tile has different label
	  label_base += layer*1000+sensor*10;
	  return label_base;
	}
      if(mms_group == mmsGroup::mms)
	{
	  // assign layer 55 and tile 0 to all
	  label_base += 55*1000 + 0;	  
	  return label_base;
	}
    }

  return -1;
}

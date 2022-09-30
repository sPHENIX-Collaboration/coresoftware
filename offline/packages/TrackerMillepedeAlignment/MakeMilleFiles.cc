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
{}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  std::string data_outfilename = "mille_output_data_file.bin";
  _mille = new Mille(data_outfilename.c_str()); 

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file("steer.txt");
  steering_file << data_outfilename.c_str() << std::endl;
  steering_file.close();

  return ret;
}

//____________________________________________________________________________..
int MakeMilleFiles::process_event(PHCompositeNode */*topNode*/)
{
  // Outline:
  //
  // loop over tracks 
  //   Make any track cuts here to skip undesirable tracks (maybe low pT?)
  //   * loop over measurements for each track
  //      * for each measurement
  //         * Get measurement value and error (global, what to use for error?)
  //         * Get track state at measurement and determine residual, error is measurement error
  //         Determine # of local parameters
  //         * Determine # of global parameters and their labels
  //         Loop over local parameters
  //            Determine derivative of residual wrt to this parameter
  //            Add to array of local derivatives, array index is used as local parameter label
  //         * Loop over global parameters
  //            * Determine derivatives of residual wrt this parameter
  //            * Add to array of global derivatives, add global par label to array of labels
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
      
      if (Verbosity() > 0)
	{
	  std::cout << std::endl << __LINE__   << ": Processing track itrack: " << phtrk_iter->first << ": nhits: " << track-> size_cluster_keys()
		    << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
	}
     
      // define running iterator over track states, used to match a given cluster to a track state
      auto state_iter = track->begin_states();

      // Get all clusters for this track
      //      std::vector<Acts::Vector3> globalClusterPositions;
      //      std::map<TrkrDefs::cluskey, Acts::Vector3> all_clusters;
      // we get the clusters from the seeds
      for( const auto& seed: { track->get_silicon_seed(), track->get_tpc_seed() } )
	{
	  for (auto key_iter = seed->begin_cluster_keys();
	       key_iter != seed->end_cluster_keys();
	       ++key_iter)
	    {
	      TrkrDefs::cluskey cluster_key = *key_iter;
	      std::cout << " looking for cluster key " << cluster_key << std::endl;
	      auto cluster = _cluster_map->findCluster(cluster_key);
	      if(!cluster) { continue;}
	      
	      // we want the global cluster positions
	      Acts::Vector3 global  = _tGeometry->getGlobalPosition(cluster_key, cluster);
	      
	      unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
	      
	      // TPC clusters need distortion corrections, silicon and MM's clusters do not
	      //=========================================================
	      if(trkrId == TrkrDefs::tpcId)
		{
		  // make all corrections to global position of TPC cluster
		  unsigned int side = TpcDefs::getSide(cluster_key);
		  unsigned int crossing = 0; // for now
		  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
		  global[2] = z;
		  
		  // apply distortion corrections
		  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
		  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
		  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
		}
	      
	      // we have our global cluster position, corrected if necessary
	      // each component of this vector is a measurement
	      // we need to find the residual and derivative in each coordinate direction
	      
	      // find track state that is the closest to cluster -  assumes that both clusters and states are sorted along r
	      float clus_radius = sqrt(global[0]*global[0]+global[1]*global[1]);
	      float dr_min = -1;
	      for( auto iter = state_iter; iter != track->end_states(); ++iter )
		{
		  const auto dr = std::abs( clus_radius - ( iter->second->get_x()* iter->second->get_x() + iter->second->get_y()* iter->second->get_y())  );
		  if( dr_min < 0 || dr < dr_min )
		    {
		      state_iter = iter;
		      dr_min = dr;
		    } else break;
		}
	      
	      // Use the PCA of the track to the measurement to get the residual
	      Acts::Vector3 pca = getPCALinePoint(global, state_iter->second);
	      Acts::Vector3 residual = global - pca;
	      
	      if(Verbosity() > 1)
		{
		  std::cout << " cluster global " << global(0) << "  " << global(1) << "  " << global(2)
			    << " track PCA " << pca(0) << "  " << pca(1) << "  " << pca(2) << std::endl;
		}
	      
	      // need standard deviation of measurements
	      float phi = atan2(global(1), global(0));
	      Acts::Vector3 clus_sigma(0,0,0);
	      if(_cluster_version==3)
		{
		  clus_sigma(2) = cluster->getZError();
		  clus_sigma(0) = cluster->getRPhiError() * cos(phi);
		  clus_sigma(1) = cluster->getRPhiError() * sin(phi);
		}
	      else if(_cluster_version==4)
		{
		  double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
		  auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster,clusRadius,cluster_key);
		  float exy2 = para_errors.first * Acts::UnitConstants::cm2;
		  float ez2 = para_errors.second * Acts::UnitConstants::cm2;
		  clus_sigma(2) = sqrt(ez2);
		  clus_sigma(0) = sqrt(exy2) * cos(phi);
		  clus_sigma(1) = sqrt(exy2) * sin(phi);
		}
	      
	      // Get the surface for this cluster
	      Surface surf = _tGeometry->maps().getSurface(cluster_key, cluster);
	      
	      // if this is a TPC cluster, check that the corrections did not change the surface
	      if(trkrId == TrkrDefs::tpcId)
		{
		  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
		  TrkrDefs::subsurfkey new_subsurfkey = 0;    
		  surf = _tGeometry->get_tpc_surface_from_coords(hitsetkey,  global, new_subsurfkey);
		}
	      if(!surf)  { continue; }

	      // For now, ignore local pars
	      static const int NLC = 0;
	      float lcl_derivative[NLC];
	      	      
	      // The global alignment parameters are given initial values of zero by default
	      // We identify the global alignment parameters for this surface
	      Acts::GeometryIdentifier id = surf->geometryId();
	      unsigned int geolayer = id.layer();
	      unsigned int sensor = id.sensitive();
	      int label_base = geolayer*1000+sensor*10;
	      
	      // Add the measurement separately for each coordinate direction to Mille
	      // set the derivatives non-zero only for parameters we want to be optimized
	      // For now, we just want to optimize the translation corrections

	      // These don't change for this track
	      static const int NGL = 6;
	      int glbl_label[NGL];
	      for(int i=0;i<NGL;++i) {glbl_label[i] = label_base + i;}
	      float glbl_derivative[NGL];

	      // x - relevant global pars are alpha, beta, gamma, dx (ipar 0,1,2,3)
	      for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
	      glbl_derivative[3] = 1.0;  // optimize dx
	      _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(0), clus_sigma(0));
	      
	      // y - relevant global pars are alpha, beta, gamma, dy (ipar 0,1,2,4)
	      for(int i=0;i<NGL;++i) {glbl_derivative[i] = 0.0;}
	      glbl_derivative[4] = 1.0; // optimize dy
	      _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(1), clus_sigma(1));
	      
	      // z - relevant global pars are alpha, beta, dz (ipar 1,1,2,5)
	      glbl_derivative[5] = 1.0;  // optimize dz
	      _mille->mille(NLC, lcl_derivative, NGL, glbl_derivative, glbl_label, residual(2), clus_sigma(2));
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
  /*
  _tpc_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_tpc_geom_container)
  {
    std::cout << PHWHERE << " ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */
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
  /*
  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( _dcc )
    { 
      std::cout << "MakeMilleFiles:   found TPC distortion correction container" << std::endl; 
    }
  */        
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 MakeMilleFiles::getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state)
{
  // Approximate track with a straight line consisting of the state position and the vector (px,py,pz)   

  Acts::Vector3 track_dir(state->get_px(), state->get_py(), state->get_pz());
  track_dir /= track_dir.norm(); 
  Acts::Vector3 track_base(state->get_x(), state->get_y(), state->get_z());

  // The position of the closest point on the line is:
  // track_base + projection of difference between the point and track_base on the line vector
  Acts::Vector3 pca = track_base + ( (global - track_base).dot(track_dir) ) * track_dir;

  return pca;
}

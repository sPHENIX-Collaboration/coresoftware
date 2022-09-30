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

//____________________________________________________________________________..
MakeMilleFiles::MakeMilleFiles(const std::string &name)
  : SubsysReco(name)
{}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output file
  std::string outfilename = "mille_output_data_file.txt";
  _mille = new Mille(outfilename.c_str()); 

  // Create global parameter labels and map of derivatives

  //MVTX & INTT
  for(auto mapIter = _tGeometry->maps().m_siliconSurfaceMap.begin(); mapIter !=  _tGeometry->maps().m_siliconSurfaceMap.end(); ++mapIter)
    {
      auto surf = mapIter->second;
      Acts::GeometryIdentifier id = surf->geometryId();
      unsigned int layer = id.layer();
      unsigned int sensor = id.sensitive();
      for(int ipar = 0; ipar < 6; ++ipar)
	{
	  int labelGL = layer * 10000 + sensor * 10 + ipar;
	  if(ipar < 3)
	    {
	      // angles
	      float derivative = 0.0;  // for now
	      derivativeGL.insert(std::make_pair(labelGL, derivative));
	    }
	  else
	    {
	      // translations
	      float derivative = 0.58;
	      derivativeGL.insert(std::make_pair(labelGL, derivative));
	    }
	}
    }

  // TPC
  for(auto mapIter = _tGeometry->maps().m_tpcSurfaceMap.begin(); mapIter !=  _tGeometry->maps().m_tpcSurfaceMap.end(); ++mapIter)
    {
      auto surf_vec = mapIter->second;
      for(auto& surf : surf_vec)
	{
	  Acts::GeometryIdentifier id = surf->geometryId();
	  unsigned int layer = id.layer();
	  unsigned int sensor = id.sensitive();
	  for(int ipar = 0; ipar < 6; ++ipar)
	    {
	      int labelGL = layer * 10000 + sensor * 10 + ipar;
	      if(ipar < 3)
		{
		  // angles
		  float derivative = 0.0;  // for now
		  derivativeGL.insert(std::make_pair(labelGL, derivative));
		}
	      else
		{
		  // translations
		  float derivative = 0.58;
		  derivativeGL.insert(std::make_pair(labelGL, derivative));
		}
	    }
	}
    }

  // TPOT
 for(auto mapIter = _tGeometry->maps().m_mmSurfaceMap.begin(); mapIter !=  _tGeometry->maps().m_mmSurfaceMap.end(); ++mapIter)
    {
      auto surf = mapIter->second;
      Acts::GeometryIdentifier id = surf->geometryId();
      unsigned int layer = id.layer();
      unsigned int sensor = id.sensitive();
      for(int ipar = 0; ipar < 6; ++ipar)
	{
	  int labelGL = layer * 10000 + sensor * 10 + ipar;
	  if(ipar < 3)
	    {
	      // angles
	      float derivative = 0.0;  // for now
	      derivativeGL.insert(std::make_pair(labelGL, derivative));
	    }
	  else
	    {
	      // translations
	      float derivative = 0.58;
	      derivativeGL.insert(std::make_pair(labelGL, derivative));
	    }
	}
    }


  // Write steering file(s) here?

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
      
      if (Verbosity() >= 1)
	{
	  std::cout << std::endl << __LINE__   << ": Processing track itrack: " << phtrk_iter->first << ": nhits: " << track-> size_cluster_keys()
		    << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
	}
      
      // running iterator over track states, used to match a given cluster to a track state
      auto state_iter = track->begin_states();

      // Get all clusters for this track
      std::vector<Acts::Vector3> globalClusterPositions;
      std::map<TrkrDefs::cluskey, Acts::Vector3> all_clusters;

      for (auto key_iter = track->begin_cluster_keys();
	   key_iter != track->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
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
	  
	  // Use the DCA of the track to the measurement as the residual
	  float dca = getDCALinePoint(global, state_iter->second);

	  // need standard deviation of measurement
	  float clus_sigma = 0.0;
	  if(_cluster_version==3){
	    clus_sigma = sqrt(cluster->getRPhiError()*cluster->getRPhiError() + cluster->getZError() * cluster->getZError());    
	  }else if(_cluster_version==4){
	    double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
	    auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster,clusRadius,cluster_key);
	    auto exy2 = para_errors.first * Acts::UnitConstants::cm2;
	    auto ez2 = para_errors.second * Acts::UnitConstants::cm2;
	    clus_sigma = sqrt(exy2+ez2);
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

	  // The global alignment parameters are given initial values of zero by default
  	  // We identify the global alignment parameters for this surface
	  Acts::GeometryIdentifier id = surf->geometryId();
	  unsigned int geolayer = id.layer();
	  unsigned int sensor = id.sensitive();

	  static const int NGL = 6;
	  int glbl_label[NGL];
	  float glbl_derivative[NGL];
	  for(int ipar=0;ipar<NGL;++ipar)
	    {
	      glbl_label[ipar] = geolayer*1000+sensor*10+ipar;
	      glbl_derivative[ipar] = derivativeGL.find(glbl_label[ipar])->second;
	    }

	  // For now
	  static const int NLC = 0;
	  float lcl_derivative[NLC];

	
	  // Add this measurement to Mille
	  _mille->mille(NLC, lcl_derivative,
			NGL, glbl_derivative, glbl_label, dca, clus_sigma);
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

float MakeMilleFiles::getDCALinePoint(Acts::Vector3 global, SvtxTrackState* state)
{
  // Approximate track with a straight line consisting of the state position and the vector (px,py,pz)   

  Acts::Vector3 track_dir(state->get_px(), state->get_py(), state->get_pz());
  Acts::Vector3 track_base(state->get_x(), state->get_y(), state->get_z());

  auto num = (global - track_base).cross(track_dir);
  float dca = num.norm() / track_dir.norm();

  return dca;
}

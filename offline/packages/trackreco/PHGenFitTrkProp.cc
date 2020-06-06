/*!
 *  \file PHGenFitTrkProp.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#include "PHGenFitTrkProp.h"

#include "AssocInfoContainer.h"

#include <trackbase/TrkrCluster.h>                      // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>                         // for cluskey, getL...

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
//

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <g4bbc/BbcVertexMap.h>

#include <phfield/PHFieldUtility.h>

#include <phgeom/PHGeomUtility.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHTimer.h>                              // for PHTimer
#include <phool/phool.h>                                // for PHWHERE

// GenFit
#include <GenFit/EventDisplay.h>                        // for EventDisplay
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>                       // for Measurement
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>

//ROOT includes for debugging
#include <TFile.h>
#include <TMatrixDSymfwd.h>                             // for TMatrixDSym
#include <TMatrixTSym.h>                                // for TMatrixTSym
#include <TMatrixTUtils.h>                              // for TMatrixTRow
#include <TNtuple.h>
#include <TVector3.h>                                   // for TVector3
#include <TVectorDfwd.h>                                // for TVectorD
#include <TVectorT.h>                                   // for TVectorT

// standard includes
#include <cassert>                                     // for assert
#include <climits>                                     // for UINT_MAX
#include <cmath>
#include <cstdlib>                                     // for exit
#include <iostream>
#include <memory>

class PHField;
class TGeoManager;
namespace genfit { class AbsTrackRep; }

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//#define _DEBUG_

//#define _USE_ALAN_FULL_VERTEXING_
#define _USE_ALAN_TRACK_REFITTING_

//#define _MEARGE_SEED_CLUSTER_
//#define _USE_ZERO_SEED_

//#define _USE_CONSTANT_SEARCH_WIN_

//#define _DO_FULL_FITTING_

using namespace std;

#ifdef _DEBUG_
ofstream fout_kalman_pull("kalman_pull.txt");
ofstream fout_chi2("chi2.txt");
#endif

PHGenFitTrkProp::PHGenFitTrkProp(
    const string& name,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc)
  : PHTrackPropagating(name)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_all(_nlayers_maps + _nlayers_intt + _nlayers_tpc)
  , _firstlayer_maps(0)
  , _firstlayer_intt(_firstlayer_maps + _nlayers_maps)
  , _firstlayer_tpc(_firstlayer_intt + _nlayers_intt)
{}

int PHGenFitTrkProp::Setup(PHCompositeNode* topNode)
{
  // Start new interface ----

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  // End new interface ----

  ret = InitializeGeometry(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  ret = InitializePHGenFit(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  if (_analyzing_mode)
  {
    cout << "Ana Mode, creating ntuples! " << endl;
    _analyzing_file = new TFile("./PHGenFitTrkProp.root", "RECREATE");
    //	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","spt:seta:sphi:pt:eta:phi:layer:ncand:nmeas");
    _analyzing_ntuple = new TNtuple("ana_nt", "ana_nt", "pt:kappa:d:phi:dzdl:z0:nhit:ml:rec:dt");
    cout << "Done" << endl;
  }

  //	ret = InitializeGeometry(topNode);
  //	if(ret != Fun4AllReturnCodes::EVENT_OK)
  //	  return ret;

  /*!
	 * Initilize parameters
	 */
  for (int layer = 0; layer < _nlayers_all; ++layer)
  {
    _search_wins_phi.insert(std::make_pair(layer, _search_win_phi));
    _search_wins_theta.insert(std::make_pair(layer, _search_win_theta));
    _max_incr_chi2s.insert(std::make_pair(layer, _max_incr_chi2));
  }

#ifdef _DEBUG_
  for (int layer = 0; layer < _nlayers_all; ++layer)
  {
    cout
        << __LINE__
        << ": layer: " << layer
        << ": search_wins_rphi: " << _search_wins_phi[layer]
        << ": search_wins_z: " << _search_wins_theta[layer]
        << ": max_incr_chi2: " << _max_incr_chi2s[layer]
        << endl;
  }
#endif

  _t_seeds_cleanup = new PHTimer("_t_seeds_cleanup");
  _t_seeds_cleanup->stop();

  _t_translate_to_PHGenFitTrack = new PHTimer("_t_translate_to_PHGenFitTrack");
  _t_translate_to_PHGenFitTrack->stop();

  _t_translate1 = new PHTimer("_t_translate1");
  _t_translate1->stop();
  _t_translate2 = new PHTimer("_t_translate2");
  _t_translate2->stop();
  _t_translate3 = new PHTimer("_t_translate3");
  _t_translate3->stop();

  _t_kalman_pat_rec = new PHTimer("_t_kalman_pat_rec");
  _t_kalman_pat_rec->stop();

  _t_search_clusters = new PHTimer("_t_search_clusters");
  _t_search_clusters->stop();

  _t_search_clusters_encoding = new PHTimer("_t_search_clusters_encoding");
  _t_search_clusters_encoding->stop();

  _t_search_clusters_map_iter = new PHTimer("_t_search_clusters_map_iter");
  _t_search_clusters_map_iter->stop();

  _t_track_propagation = new PHTimer("_t_track_propergation");
  _t_track_propagation->stop();

  _t_full_fitting = new PHTimer("_t_full_fitting");
  _t_full_fitting->stop();

  _t_output_io = new PHTimer("_t_output_io");
  _t_output_io->stop();

  return ret;
}

int PHGenFitTrkProp::Process()
{
  if (Verbosity() > 10)
  {
    cout << "PHGenFitTrkProp::process_event -- entered" << endl;
    cout << "nMapsLayers = " << _nlayers_maps << endl;
    cout << "nInttLayers = " << _nlayers_intt << endl;
    cout << "nTPCLayers = " << _nlayers_tpc << endl;
  }
  // start fresh
  _gftrk_hitkey_map.clear();

  _vertex.clear();
  _vertex_error.clear();

  if (_vertex_map)
    {
      for(unsigned int ivert=0; ivert<_vertex_map->size(); ++ivert)  
	{
	  SvtxVertex* vertex = _vertex_map->get(ivert);

	  std::vector<float> v = { vertex->get_x(), vertex->get_y(), vertex->get_z() };
	  _vertex.push_back(v);

	  std::vector<float> v_err = { sqrt(vertex->get_error(0,0)), sqrt(vertex->get_error(1,1)), sqrt(vertex->get_error(2,2)) }; 
	  _vertex_error.push_back(v_err);
	}
    }
  
  {
    //-----------------------------------
    // Kalman track propagating
    //-----------------------------------
    if (Verbosity() >= 1) _t_kalman_pat_rec->restart();
    int ret = KalmanTrkProp();
    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
    if (Verbosity() >= 1) _t_kalman_pat_rec->stop();
  }
  if (Verbosity() > 1) print_timers();

  ++_event;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGenFitTrkProp::print_timers()
{
  std::cout << "=============== PHGenFitTrkProp::print_timers: ===============" << std::endl;
  std::cout << "\t - Seeds Cleanup:          " << _t_seeds_cleanup->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "CPUSCALE Pattern recognition time:    " << _t_kalman_pat_rec->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Track Translation time: " << _t_translate_to_PHGenFitTrack->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation1 time: " << _t_translate1->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation2 time: " << _t_translate2->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation3 time: " << _t_translate3->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Cluster searching time: " << _t_search_clusters->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t\t - Encoding time:        " << _t_search_clusters_encoding->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t\t - Map iteration:        " << _t_search_clusters_map_iter->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Kalman updater time:    " << _t_track_propagation->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "Full fitting time:           " << _t_full_fitting->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "Output IO time:              " << _t_output_io->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "=======================================" << std::endl;
}

int PHGenFitTrkProp::End()
{
  if (_do_evt_display)
    _fitter->displayEvent();

  delete _t_seeds_cleanup;
  delete _t_translate_to_PHGenFitTrack;
  delete _t_translate1;
  delete _t_translate2;
  delete _t_translate3;
  delete _t_kalman_pat_rec;
  delete _t_full_fitting;
  delete _t_search_clusters;
  delete _t_track_propagation;
  delete _t_output_io;

#ifdef _DEBUG_
  LogDebug("Leaving End \n");
#endif

#ifdef _DEBUG_
  fout_kalman_pull.close();
  fout_chi2.close();
#endif

  if (_analyzing_mode)
  {
    cout << " cleaning up " << endl;
    _analyzing_file->cd();
    _analyzing_ntuple->Write();
    _analyzing_file->Close();
    //	  delete _analyzing_ntuple;
    // delete _analyzing_file;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::InitializePHGenFit(PHCompositeNode* topNode)
{
  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);

  PHField* field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

  _fitter.reset( PHGenFit::Fitter::getInstance(tgeo_manager, field, _track_fitting_alg_name, "RKTrackRep", _do_evt_display) ); 

  if (!_fitter)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

#ifdef _DEBUG_
  _fitter->set_verbosity(10);
#endif

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  _topNode = topNode;

  // used in fast vertexing from BBC
  _bbc_vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");

  // _cluster_map is the TrkrClusterMap container, and is found in PHTrackPropagating
 
  _geom_container_intt = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  _geom_container_maps = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::check_track_exists(MapPHGenFitTrack::iterator iter, SvtxTrackMap::Iter phtrk_iter)
{
  //Loop over hitIDs on current track and check if they have been used
  unsigned int n_clu = iter->second->get_cluster_keys().size();
  
  unsigned int n_clu_used = 0;

  const std::vector<TrkrDefs::cluskey>& clusterkeys = iter->second->get_cluster_keys();

  int n = 0;
  for (TrkrDefs::cluskey iCluId = 0; iCluId < clusterkeys.size(); ++iCluId)
  {
    TrkrDefs::cluskey cluster_ID = clusterkeys[iCluId];
    
    if(_gftrk_hitkey_map.count(iCluId)>0) n_clu_used++;
    if (Verbosity() >= 10){
      cout << " trk map size: " << _gftrk_hitkey_map.count(iCluId) << endl;
      cout << "n: " << n << "#Clu_g = " << iCluId << " ntrack match: "  << _assoc_container->GetTracksFromCluster(cluster_ID).size()
	   << " layer: " << (float)TrkrDefs::getLayer(cluster_ID)
	   << " r: " << sqrt(_cluster_map->findCluster(cluster_ID)->getX()*_cluster_map->findCluster(cluster_ID)->getX() +_cluster_map->findCluster(cluster_ID)->getY()*_cluster_map->findCluster(cluster_ID)->getY() )
	   << " used: " << n_clu_used
	   << endl;
      n++;
    }
  }

  int code = 0;
  if (((float) n_clu_used / n_clu) > 0.3)
  {
    if (Verbosity() >= 1)
      cout << "Found duplicate track. n_clu: " << n_clu << " c_clu_used: " << n_clu_used << endl;
    /*
    for(TrkrDefs::cluskey iCluId = 0; iCluId < clusterkeys.size(); ++iCluId){
      TrkrDefs::cluskey cluster_ID = clusterkeys[iCluId];
     
    }
    */
    
    return code;
  }
  code = 1;
  return code;
}

int PHGenFitTrkProp::KalmanTrkProp()
{
  // _track_map is created in PHHoughSeeding::ExportOutput at line 1714 and written to SvtxTrack node.
  // PHTrackPropagating (base class of this one) gets it from the node tree and calls this "Process" method.
  // The process method then calls this method.

#ifdef _DEBUG_
  std::cout << "=========================" << std::endl;
  std::cout << "PHGenFitTrkProp::KalmanTrkProp: Start: Event: " << _event << std::endl;
  std::cout << "Total Raw Tracks: " << _track_map->size() << std::endl;
  std::cout << "=========================" << std::endl;
#endif

  /*!
   *   sort clusters
   */  
  // A map is needed for each vertex location. This is a vector of maps
  _layer_thetaID_phiID_clusterID.clear();
  /*  for(unsigned int ivert=0; ivert<_vertex.size(); ++ivert)
    {
      BuildLayerZPhiHitMap(ivert);      
    }
  */
  BuildLayerZPhiHitMap(0);

  vector<genfit::Track*> evt_disp_copy;

  _PHGenFitTracks.clear();

  //_track_map->identify();

  if (Verbosity() > 1){
    cout << " found " << _track_map->size() << " track seeds " << endl;
  }

  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end();)
    {
      SvtxTrack* tracklet = phtrk_iter->second;
      if (Verbosity() >= 1){
	std::cout
	  << __LINE__
	  << ": Processing seed itrack: " << phtrk_iter->first
	  << ": nhits: " << tracklet-> size_cluster_keys()
	  << ": Total tracks: " << _track_map->size()
	  << endl;
      }
      _ntrack = phtrk_iter->first;
      /*!
       * Translate sPHENIX track To PHGenFitTracks
       */
      if (Verbosity() > 1) _t_translate_to_PHGenFitTrack->restart();
      /*
      std::shared_ptr<PHGenFit::Track> rf_phgf_track =  ReFitTrack(tracklet);
      if(!rf_phgf_track)
	{
	  cout << "FAIL ReFitting input track# " << phtrk_iter->first << endl;
	}
      */
      //if (Verbosity() > 1) _t_translate_to_PHGenFitTrack->stop();
      
      /*!
       * Handle track propagation, termination, output and evt disp.
       */
      bool is_splitting_track = false;
#ifdef _DEBUG_
      int i = 0;
#endif
      SvtxTrackToPHGenFitTracks(tracklet);
      if (_PHGenFitTracks.empty()) 
	{
	  cout << "FAIL SvtxTrackToGenFit " <<  phtrk_iter->first << endl;
	  //cout << "Warning: Conversion of SvtxTrack tracklet " <<  phtrk_iter->first << " to PHGenFitTrack failed, moving to next tracklet " << endl;
	  ++phtrk_iter;
	  continue;
	}
      
      for (auto gftrk_iter = _PHGenFitTracks.begin();
	   gftrk_iter != _PHGenFitTracks.end(); ++gftrk_iter)
	{
	  if(Verbosity() > 10)
	    {
	      cout
		<< __LINE__
		<< ": propagating Genfit track: " << phtrk_iter->first << endl;
	    }

	  // associate this track with the same vertex as the seed track
	  //	  unsigned int ivert = tracklet->get_vertex_id();
	  unsigned int ivert = 0;
	  gftrk_iter->second->set_vertex_id(ivert);
	  //cout << PHWHERE << " read back track vertex ID from Genfit track " << gftrk_iter->second->get_vertex_id() << endl;
	  if(ivert > _vertex.size())
	    {
	      cout << PHWHERE << " Track vertex is screwed up, have to quit! #" << ivert << " v_size: " << _vertex.size()  << endl;
	      return Fun4AllReturnCodes::ABORTRUN;
	    }
	  
	  std::vector<TrkrDefs::cluskey> clusterkeys = gftrk_iter->second->get_cluster_keys();
	  
	  unsigned int init_layer = UINT_MAX;
	  
	  if (!is_splitting_track)
	    {
	      if (_init_direction == 1)
		{
		  //init_layer = _cluster_map->get(clusterkeys.front())->get_layer();
		  init_layer = TrkrDefs::getLayer(clusterkeys.front());
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, true);
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, false);
		}
	      else
		{
		  //init_layer = _cluster_map->get(clusterkeys.back())->get_layer();
		  init_layer = TrkrDefs::getLayer(clusterkeys.back());
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, true);
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, false);
		}
	      is_splitting_track = true;
	    }
	  else
	    {
	      if (_init_direction == 1)
		{
		  //init_layer = _cluster_map->get(clusterkeys.front())->get_layer();
		  init_layer = TrkrDefs::getLayer(clusterkeys.front());
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, false);
		}
	      else
		{
		  //init_layer = _cluster_map->get(clusterkeys.back())->get_layer();
		  init_layer = TrkrDefs::getLayer(clusterkeys.back());
		  TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, false);
		}
	    }
	  
#ifdef _DEBUG_
	  cout
	    << __LINE__
	    << ": tracki: " << i
	    << ": clusterkeys size:  " << gftrk_iter->second->get_cluster_keys().size()
	    << ": quality: " << gftrk_iter->first
	    << endl;
	  ++i;
#endif
	  
	  //_trackID_PHGenFitTrack.erase(iter);
	}  // loop _PHGenFitTracks
      
      
#ifdef _DEBUG_
      i = 0;
      for (auto iter = _PHGenFitTracks.begin();
	   iter != _PHGenFitTracks.end(); ++iter)
	{
	  cout
	    << __LINE__
	    << ": track: " << i++
	    << ": clusterkeys size:  " << iter->second->get_cluster_keys().size()
	    << ": quality: " << iter->first
	    << endl;
	}
#endif
      
      //std::sort(_PHGenFitTracks.begin(), _PHGenFitTracks.end());
      _PHGenFitTracks.sort();
      
#ifdef _DEBUG_
      for (auto iter = _PHGenFitTracks.begin();
	   iter != _PHGenFitTracks.end(); ++iter)
	{
	  cout
	    << __LINE__
	    << ": clusterkeys size:  " << iter->second->get_cluster_keys().size()
	    << ": quality: " << iter->first
	    << endl;
	}
#endif
      
      auto gftrk_iter_best = _PHGenFitTracks.begin();
      
      int track_exists = check_track_exists(gftrk_iter_best,phtrk_iter);
      
      if (gftrk_iter_best->second->get_cluster_keys().size() >= _min_good_track_hits && track_exists)
	{
	  OutputPHGenFitTrack(gftrk_iter_best, phtrk_iter);
#ifdef _DEBUG_
	  cout << __LINE__ << endl;
#endif
	  if (_do_evt_display)
	    {
	      evt_disp_copy.push_back(
				      new genfit::Track(*gftrk_iter_best->second->getGenFitTrack()));
	    }
	}
      else
	{
	  auto key = phtrk_iter->first;
	  ++phtrk_iter;
	  _track_map->erase(key);
	  continue;
	}
      
      _PHGenFitTracks.clear();
      
      ++phtrk_iter;
    }
  
  if(Verbosity() > 1)
    {
      std::cout << "=========================" << std::endl;
      std::cout << "PHGenFitTrkProp::KalmanTrkProp: End: Event: " << _event << std::endl;
      std::cout << "PHGenFitTrkProp::KalmanTrkProp: End: Event: " << _event << std::endl;
      std::cout << "Total Final Tracks: " << _track_map->size() << std::endl;
    }
  
  if (_do_evt_display)
    {
      _fitter->getEventDisplay()->addEvent(evt_disp_copy);
    }
  else
    {
      evt_disp_copy.clear();
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::OutputPHGenFitTrack(
					 MapPHGenFitTrack::iterator gftrk_iter,
					 SvtxTrackMap::Iter phtrk_iter)
{
  if(Verbosity() > 10)
    {
      std::cout << "=========================" << std::endl;
      std::cout << __LINE__ << ": iPHGenFitTrack: " << phtrk_iter->first << std::endl;
      std::cout << __LINE__ << ": _track_map->size(): " << _track_map->size() << std::endl;
      std::cout << "Contains: " << gftrk_iter->second->get_cluster_keys().size() << " clusters." << std::endl;
      std::cout << "=========================" << std::endl;
    }

  // get the track from the node tree and extract the track and vertex id  
  SvtxTrack* track = phtrk_iter->second;
  auto track_id = track->get_id();
  auto vertex_id = track->get_vertex_id();

  // now reset the track info and rewrite it using info from the Genfit track
  track->Reset();
  track->set_id(track_id);
  track->set_vertex_id(vertex_id);
  
#ifdef _DO_FULL_FITTING_
  if (Verbosity() >= 1) _t_full_fitting->restart();
  if (_fitter->processTrack(gftrk_iter->second.get(), false) != 0)
    {
      if (Verbosity() >= 1)
	LogWarning("Track fitting failed\n");
      //delete track;
      return -1;
    }
  if (Verbosity() >= 1) _t_full_fitting->stop();
  
  if (Verbosity() >= 1) _t_output_io->restart();
  //		iter->second->getGenFitTrack()->Print();
  
  track.set_chisq(gftrk_iter->second->get_chi2());
  track.set_ndf(gftrk_iter->second->get_ndf());
  
  // Use fitted vertex
  TVector3 vertex_position(_vertex[vertex_id][0], _vertex[vertex_id][1], _vertex[vertex_id][2]);
  std::unique_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca = nullptr;
  try
    {
      gf_state_vertex_ca = std::unique_ptr<genfit::MeasuredStateOnPlane>(gftrk_iter->second->extrapolateToPoint(vertex_position));
    }
  catch (...)
    {
      if (Verbosity() >= 2)
	LogWarning("extrapolateToPoint failed!");
    }
  if (!gf_state_vertex_ca)
    {
      //delete out_track;
      return -1;
    }
  
  TVector3 mom = gf_state_vertex_ca->getMom();
  TVector3 pos = gf_state_vertex_ca->getPos();
  TMatrixDSym cov = gf_state_vertex_ca->get6DCov();

#else
  TVectorD state = gftrk_iter->second->getGenFitTrack()->getStateSeed();
  TVector3 pos(state(0), state(1), state(2));
  TVector3 mom(state(3), state(4), state(5));
  TMatrixDSym cov = gftrk_iter->second->getGenFitTrack()->getCovSeed();
#endif

  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      track->set_error(i, j, cov[i][j]);
    }
  }

  track->set_px(mom.Px());
  track->set_py(mom.Py());
  track->set_pz(mom.Pz());
  
  track->set_x(pos.X());
  track->set_y(pos.Y());
  track->set_z(pos.Z());
  
  for (TrkrDefs::cluskey cluster_key : gftrk_iter->second->get_cluster_keys())
    {
      if(Verbosity() > 10) cout << " track id: " << phtrk_iter->first <<  " adding clusterkey " << cluster_key << endl;
      _gftrk_hitkey_map.insert(std::make_pair(cluster_key, phtrk_iter->first));
      track->insert_cluster_key(cluster_key);
    }
  
  //Check track quality
  //		bool is_good_track = true;
  
  Int_t n_maps = 0;
  Int_t n_intt = 0;
  Int_t n_tpc = 0;
  
  for (auto iter = track->begin_cluster_keys(); iter != track->end_cluster_keys(); ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    unsigned int layer = TrkrDefs::getLayer(cluster_key);
    if( is_maps_layer( layer ) )
    {
      n_maps++;
    } else if( is_intt_layer( layer ) ) {
      n_intt++;
      if (n_intt > 4)
      {
        cout << PHWHERE << " Can not have more than 4 INTT layers, quit!" << endl;
        exit(1);
      }
    } else if( is_tpc_layer( layer ) ) {
      n_tpc++;
    }
  }
  
  // Add the cluster-track association to the association table for later use
  for (TrkrDefs::cluskey cluster_key : gftrk_iter->second->get_cluster_keys())
    {
      //cout << "Add cluster key " << cluster_key << " to ClusterTrackAssoc " << track->get_id() << endl;
      _assoc_container->SetClusterTrackAssoc(cluster_key, track->get_id());
    }

  if (Verbosity() > 10)
    {
      cout << "Output propagared track " << track->get_id() << " vertex " << track->get_vertex_id()<< " quality = " << track->get_quality()
	   << " clusters: mvtx " << n_maps << " intt " << n_intt << " tpc  " << n_tpc 
	   << endl;
      cout << "px = " << track->get_px() << " py = " << track->get_py()
	   << " pz = " << track->get_pz() << endl;
    }
  
  if (Verbosity() >= 1) _t_output_io->stop();
  
  return 0;
}

std::shared_ptr<PHGenFit::Track> PHGenFitTrkProp::ReFitTrack(SvtxTrack* intrack)
{
  if (!intrack)
  {
    cerr << PHWHERE << " Input SvtxTrack is nullptr!" << endl;
    return nullptr;
  }
  
  PHG4CylinderGeomContainer* geom_container_intt = findNode::getClass<
      PHG4CylinderGeomContainer>(_topNode, "CYLINDERGEOM_INTT");

  PHG4CylinderGeomContainer* geom_container_mvtx = findNode::getClass<
      PHG4CylinderGeomContainer>(_topNode, "CYLINDERGEOM_MVTX");

  // prepare seed
  TVector3 seed_mom(100, 0, 0);
  TVector3 seed_pos(0, 0, 0);
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
    }
  }

  //Sort hits in r
  std::map<float, TrkrDefs::cluskey> m_r_cluster_id;
  for (auto iter = intrack->begin_cluster_keys();
       iter != intrack->end_cluster_keys(); ++iter){
    TrkrDefs::cluskey cluster_key = *iter;
    TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
    float x = cluster->getPosition(0);
    float y = cluster->getPosition(1);
    float r = sqrt(x*x + y*y);
    m_r_cluster_id.insert(std::pair<float, TrkrDefs::cluskey>(r, cluster_key));
    int layer_out = TrkrDefs::getLayer(cluster_key);
    if(Verbosity() >= 2) cout << "    Layer " << layer_out << " cluster " << cluster_key << " radius " << r << endl;
  }

  // Create measurements
  std::vector<PHGenFit::Measurement*> measurements;
  for (auto iter = m_r_cluster_id.begin();
       iter != m_r_cluster_id.end();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->second;
    const int layer = TrkrDefs::getLayer(cluster_key);

    TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
    if (!cluster){
      LogError("No cluster Found!");
      continue;
    }
    //    PHGenFit::Measurement* meas = TrkrClusterToPHGenFitMeasurement(cluster);
    
    TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

    seed_mom.SetPhi(pos.Phi());
    seed_mom.SetTheta(pos.Theta());

    //TODO use u, v explicitly?
    TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);
    
    // get the trkrid
    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);

    if(trkrid == TrkrDefs::mvtxId)
      {
	int stave_index = MvtxDefs::getStaveId(cluster_key);
	int chip_index = MvtxDefs::getChipId(cluster_key);
	
	double ladder_location[3] = {0.0, 0.0, 0.0};
	CylinderGeom_Mvtx* geom =
          dynamic_cast<CylinderGeom_Mvtx*>(geom_container_mvtx->GetLayerGeom(layer));
	// returns the center of the sensor in world coordinates - used to get the ladder phi location
	geom->find_sensor_center(stave_index, 0,
				 0, chip_index, ladder_location);
	
	//cout << " MVTX stave phi tilt = " <<  geom->get_stave_phi_tilt()
	//   << " seg.X " << ladder_location[0] << " seg.Y " << ladder_location[1] << " seg.Z " << ladder_location[2] << endl;
	n.SetXYZ(ladder_location[0], ladder_location[1], 0);
	n.RotateZ(geom->get_stave_phi_tilt());
      }
    else if(trkrid == TrkrDefs::inttId)
      {
	CylinderGeomIntt* geom =
          dynamic_cast<CylinderGeomIntt*>(geom_container_intt->GetLayerGeom(layer));
	double hit_location[3] = {0.0, 0.0, 0.0};
	geom->find_segment_center(InttDefs::getLadderZId(cluster_key),
				  InttDefs::getLadderPhiId(cluster_key), hit_location);
	
	//cout << " Intt strip phi tilt = " <<  geom->get_strip_phi_tilt()
	//   << " seg.X " << hit_location[0] << " seg.Y " << hit_location[1] << " seg.Z " << hit_location[2] << endl;
	n.SetXYZ(hit_location[0], hit_location[1], 0);
	n.RotateZ(geom->get_strip_phi_tilt());
      }
    // end new
    //-----------------
    
    PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
								  cluster->getRPhiError(), cluster->getZError());
    cout << "Add meas layer " << layer << " cluskey " << cluster_key
	 << endl
	 << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
	 << " RPhiErr " << cluster->getRPhiError()
	 << " ZErr " << cluster->getZError()
	 << endl;
    /*
    if(Verbosity() > 10)
      {
	cout << "Add meas layer " << layer << " cluskey " << cluster_key
	     << endl
	     << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
	     << "  n.X " <<  n.X() << " n.Y " << n.Y()
	     << " RPhiErr " << cluster->getRPhiError()
	     << " ZErr " << cluster->getZError()
	     << endl;
      }
    */
    measurements.push_back(meas);
  }
  
  /*!
   * mu+:	-13
   * mu-:	13
   * pi+:	211
   * pi-:	-211
   * e-:	11
   * e+:	-11
   */
  //TODO Add multiple TrackRep choices.
  //int pid = 211;

  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
  std::shared_ptr<PHGenFit::Track> track(new PHGenFit::Track(rep, seed_pos, seed_mom,
                                                             seed_cov));

  //TODO unsorted measurements, should use sorted ones?
  track->addMeasurements(measurements);
  if (Verbosity() >= 1)
    cout <<  "#trk: " << _ntrack << " ReFit nhits: " <<  measurements.size() << endl;
  /*!
   *  Fit the track
   *  ret code 0 means 0 error or good status
   */
  
  if (_fitter->processTrack(track.get(), false) != 0)
    {
      if (Verbosity() >= 1)
	{
	  cout <<  "#trk: " << _ntrack << "Track fitting failed (ReFit) nhits: " <<  measurements.size() << endl;
	  //	  LogWarning("Track fitting failed ");
	  cout << " " << endl;
	  /*
	    cout << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
	    << " mom.X " << track->get_mom().X()
	    << " mom.Y " << track->get_mom().Y()
	    << " mom.Z " << track->get_mom().Z()
	    << endl;
	  */
	}
      //delete track;
      return nullptr;
    }
  return track;
}

int PHGenFitTrkProp::SvtxTrackToPHGenFitTracks(const SvtxTrack* svtxtrack)
{
  // clean up working array for each event
  _PHGenFitTracks.clear();

  TVector3 seed_mom(1000, 0, 0);
  TVector3 seed_pos(0, 0, 0);
  //  const float blowup_factor = 1.;
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
      //      seed_cov[i][j] = blowup_factor * svtxtrack->get_error(i, j);
    }
  }


  //Sort hits in r 
  //  std::multimap<float, TrkrDefs::cluskey> m_r_clusterID;
  std::map<float, TrkrDefs::cluskey> m_r_clusterID;
  for (auto hit_iter = svtxtrack->begin_cluster_keys();
       hit_iter != svtxtrack->end_cluster_keys(); ++hit_iter)
  {
    TrkrDefs::cluskey clusterkey = *hit_iter;    
    TrkrCluster *cluster = _cluster_map->findCluster(clusterkey);
    float r = sqrt(
		   cluster->getPosition(0) * cluster->getPosition(0) +
		   cluster->getPosition(1) * cluster->getPosition(1));
    m_r_clusterID.insert(std::pair<float, TrkrDefs::cluskey>(r, clusterkey)); 
    if (Verbosity() >= 10){
      cout << PHWHERE << " inserted r " << r << " clusterkey " << clusterkey 
	   << " layer: " << (float)TrkrDefs::getLayer(clusterkey)
	   << " max: " << svtxtrack->size_cluster_keys()
	   <<  endl;
    }
  }

  // Create measurements
  std::vector<PHGenFit::Measurement*> measurements;
  for (auto iter = m_r_clusterID.begin();
       iter != m_r_clusterID.end();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->second;

    TrkrCluster *cluster = _cluster_map->findCluster(cluster_key);
    if (!cluster){
      LogError("No cluster Found!\n");
      continue;
    }
    //    const int layer = TrkrDefs::getLayer(cluster_key);
    PHGenFit::Measurement* meas = TrkrClusterToPHGenFitMeasurement(cluster);
    TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));
    seed_mom.SetPhi(pos.Phi());
    seed_mom.SetTheta(pos.Theta());
    /*
    cout << "Add meas layer " << layer << " cluskey " << cluster_key
	 << endl
	 << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
	 << " RPhiErr " << cluster->getRPhiError()
	 << " ZErr " << cluster->getZError()
	 << endl;
    */
    if (meas){
      measurements.push_back(meas);
    }
  }

  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
  std::shared_ptr<PHGenFit::Track> track(
      new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov));

  track->addMeasurements(measurements);
  if (Verbosity() >= 1)
    cout << "#trk: " << _ntrack << " ToPHGen size:" << measurements.size()  <<endl;

  if (_fitter->processTrack(track.get(), false) != 0)
  {
    if (Verbosity() >= 1){
      cout << "#trk: " << _ntrack << " Seed fitting failed (ToPHGen) size:" << measurements.size()  <<endl;
    }
    return -1;
  }

  int nhits = track->get_cluster_keys().size();
  float chi2 = track->get_chi2();
  float ndf = track->get_ndf();

  if (nhits > 0 and chi2 > 0 and ndf > 0)
  {
    _PHGenFitTracks.push_back(
        MapPHGenFitTrack::value_type(
            PHGenFitTrkProp::TrackQuality(nhits, chi2, ndf, nhits, 0, 0), track));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::TrackPropPatRec(
				     const unsigned int ivert,
				     MapPHGenFitTrack::iterator& track_iter,
				     unsigned int init_layer, unsigned int end_layer,
				     const bool use_fitted_state_once)
{
#ifdef _DEBUG_
  cout
      << __LINE__
      << " TrackPropPatRec"
      << " : init_layer: " << init_layer
      << " : end_layer: " << end_layer
      << " : use_fitted_state_once: " << use_fitted_state_once
      << endl;
  #endif

  std::shared_ptr<PHGenFit::Track>& track = track_iter->second;
  if(Verbosity() > 10) cout << endl << PHWHERE << " Entering TrackPropPatRec for track with vertex ID " << ivert << endl;
  if(ivert >= _vertex.size())
    {
      cout << PHWHERE << " WARNING: vertex ID out of range, something wrong, quit! " << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  int direction = end_layer >= init_layer ? 1 : -1;
  assert(direction == 1 or direction == -1);

  int first_extrapolate_base_TP_id = -1;

  bool use_fitted_state = use_fitted_state_once;
  float blowup_factor = use_fitted_state ? _blowup_factor : 1.;
  unsigned int layer_occupied[_nlayers_all];
  for(int i = 0;i<_nlayers_all;i++) layer_occupied[i] = 0;
  /*!
	 * Find the last layer of with TrackPoint (TP)
	 * Asumming measuremnts are sorted by radius
	 * and cluster keys are syncronized with the TP IDs
	 */
  
  std::vector<TrkrDefs::cluskey> clusterkeys = track->get_cluster_keys();
  
  for (unsigned int i = 0; i < clusterkeys.size(); ++i)
    {
      unsigned int layer = TrkrDefs::getLayer(clusterkeys[i]);
      layer_occupied[layer] = 1;
      if(layer == init_layer)
      {
        first_extrapolate_base_TP_id = i;
        break;
      }
    }


  if (first_extrapolate_base_TP_id < 0)
  {
    if (Verbosity() > 0)
      LogError("first_extrapolate_base_TP_id < 0");
    return -1;
  }

  int extrapolate_base_TP_id = first_extrapolate_base_TP_id;

  unsigned int consecutive_missing_layer = 0;

  //	unsigned int layer = init_layer + direction;
  //	while (layer>=0 and layer < (unsigned int)_nlayers_all and layer!=end_layer) {

  for (unsigned int layer = init_layer + direction;
       layer != end_layer + direction;
       layer += direction)
  {
    // layer is unsigned int, check for >=0 is meaningless
    if (layer >= (unsigned int) _nlayers_all) break;
    if (layer_occupied[layer]) continue;

    /*!
		 * if miss too many layers terminate track propagating
		 */
    if (consecutive_missing_layer > _max_consecutive_missing_layer)
    {
      if (Verbosity() > 1)
      {
        LogWarning("consecutive_missing_layer > ") << _max_consecutive_missing_layer << endl;
      }
      if (track->get_cluster_keys().size() >= _min_good_track_hits)
        return 0;
      else
        return -1;
    }

    bool layer_updated = false;

    float layer_r = _radii_all[_layer_ilayer_map_all[layer]];

    if(Verbosity() > 10)
      {
	std::cout << "=========================" << std::endl;
	std::cout << __LINE__ << ": Event: " << _event << ": _PHGenFitTracks.size(): " << _PHGenFitTracks.size() << ": layer: " << layer << std::endl;
	std::cout << "=========================" << std::endl;
      }

#ifdef _DEBUG_
    {
      unsigned int tempIdx = extrapolate_base_TP_id >= 0 ? extrapolate_base_TP_id : extrapolate_base_TP_id + track->get_cluster_keys().size();
      cout
          << __LINE__
          << " tempIdx: " << tempIdx
          << endl;
      // tempIdx is unsigned int, checking for >=0 is meaningless
      if (tempIdx < track->get_cluster_keys().size())
      {
	TrkrDefs::cluskey extrapolate_base_cluster_id = track->get_cluster_keys()[tempIdx];
        //SvtxCluster* extrapolate_base_cluster = _cluster_map->get(extrapolate_base_cluster_id);
	TrkrCluster* extrapolate_base_cluster = _cluster_map->findCluster(extrapolate_base_cluster_id);
	int from_layer = TrkrDefs::getLayer(extrapolate_base_cluster_id);
        cout
            << __LINE__
            << ": Target layer: { " << layer
            << ", " << layer_r
            << "} : From layer: { " << from_layer
            << ", " << sqrt(extrapolate_base_cluster->getX() * extrapolate_base_cluster->getX() + extrapolate_base_cluster->getY() * extrapolate_base_cluster->getY())
            << "} : ID: " << extrapolate_base_cluster_id
            << endl;
      }
    }
#endif

    //		bool have_tp_with_fit_info = false;
    //		std::vector<unsigned int> clusterkeys = track->get_cluster_keys();
    //		for (unsigned int i = clusterkeys.size() - 1; i >= 0; --i) {
    //			std::unique_ptr<genfit::MeasuredStateOnPlane> kfsop = nullptr;
    //			genfit::Track genfit_track = track->getGenFitTrack();
    //			if (genfit_track->getNumPointsWithMeasurement() > 0) {
    //
    //				genfit::TrackPoint* tp = genfit_track->getPointWithMeasurement(tr_point_id);
    //				if (dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))) {
    //
    //					if (static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getForwardUpdate()) {
    //						have_tp_with_fit_info = true;
    //						break;
    //					}
    //				}
    //			}
    //		}

    std::unique_ptr<genfit::MeasuredStateOnPlane> state = nullptr;
    try
    {
      state = std::unique_ptr<genfit::MeasuredStateOnPlane>(track->extrapolateToCylinder(layer_r, TVector3(0, 0, 0),
                                                                                         TVector3(0, 0, 1), extrapolate_base_TP_id, direction));
      //		genfit::MeasuredStateOnPlane *state = track->extrapolateToCylinder(
      //				layer_r, TVector3(0, 0, 0), TVector3(0, 0, 1), 0);
    }
    catch (...)
    {
      if (Verbosity() > 1)
      {
        LogWarning("Can not extrapolate to Cylinder!") << std::endl;
      }
      continue;
    }

    if (!state)
    {
      if (Verbosity() > 1)
      {
        LogWarning("Can not extrapolate to Cylinder!") << std::endl;
      }
      continue;
    }

    // The vertex position is used here to caculate theta and phi for the
    // extrapolated track position. These values of theta and phi are compared 
    // with the cluster theta and phi positions in the map  "_layer_thetaID_phiID_clusterID"

    // Need to use the correct vertex position for this track as well as for the theta-phi map.

    TVector3 pos = state->getPos();
    if(Verbosity() > 10) cout << PHWHERE << " pos z " << pos.Z() << " _vertex z " << _vertex[ivert][2] << endl;
    pos.SetXYZ(
        pos.X() - _vertex[ivert][0],
        pos.Y() - _vertex[ivert][1],
        pos.Z() - _vertex[ivert][2]);

    float phi_center = pos.Phi();
    float theta_center = pos.Theta();

#ifdef _USE_CONSTANT_SEARCH_WIN_

    float phi_window = 25e-4;
    float theta_window = 25e-4;

    if (layer >= 3 and layer <= 6)
    {
      phi_window = 300e-4;
      theta_window = 0.2;
    }

    if (layer <= 2)
    {
      phi_window = 3000e-4;
      theta_window = 3000e-4;
    }
#else
    TMatrixDSym cov = state->get6DCov();

    float phi_window = _search_wins_phi[layer] * sqrt(cov[0][0] + cov[1][1] + cov[0][1] + cov[1][0]) / pos.Perp();
    float theta_window = _search_wins_theta[layer] * sqrt(cov[2][2]) / pos.Perp();

    if( is_maps_layer( layer ) )
    {
      if (phi_window > _max_search_win_phi_maps) phi_window = _max_search_win_phi_maps;
      if (phi_window < _min_search_win_phi_maps) phi_window = _min_search_win_phi_maps;
      if (theta_window > _max_search_win_theta_maps) theta_window = _max_search_win_theta_maps;
      if (theta_window < _min_search_win_theta_maps) theta_window = _min_search_win_theta_maps;
    }
    else if( is_intt_layer( layer ) )
    {
      const auto layer_intt = layer - _firstlayer_intt;
      if (phi_window > _max_search_win_phi_intt[layer_intt]) phi_window = _max_search_win_phi_intt[layer_intt];
      if (phi_window < _min_search_win_phi_intt[layer_intt]) phi_window = _min_search_win_phi_intt[layer_intt];
      if (theta_window > _max_search_win_theta_intt[layer_intt]) theta_window = _max_search_win_theta_intt[layer_intt];
      if (theta_window < _min_search_win_theta_intt[layer_intt]) theta_window = _min_search_win_theta_intt[layer_intt];
    }
    else if( is_tpc_layer( layer ) )
    {
      if (phi_window > _max_search_win_phi_tpc) phi_window = _max_search_win_phi_tpc;
      if (phi_window < _min_search_win_phi_tpc) phi_window = _min_search_win_phi_tpc;
      if (theta_window > _max_search_win_theta_tpc) theta_window = _max_search_win_theta_tpc;
      if (theta_window < _min_search_win_theta_tpc) theta_window = _min_search_win_theta_tpc;
    }

    //FIXME optimize this
    //		if(layer == _nlayers_maps + _nlayers_intt -1) {
    //			phi_window = 0.02;
    //			theta_window = 0.04;
    //		}

#endif

#ifdef _DEBUG_
    cout << __LINE__ << ": ";
    printf("layer: %u: r: %f: phi: %f +- %f; theta: %f +- %f\n",
           layer, pos.Perp(),
           phi_center, phi_window,
           theta_center, theta_window);

//		cout<<__LINE__<<": ";
//		printf("layer: %d:  phi: %f +- %f\n",
//				layer,
//				pos.Phi(), phi_window
//				);
#endif

    if (Verbosity() >= 1) _t_search_clusters->restart();
    std::vector<TrkrDefs::cluskey> new_cluster_keys = SearchHitsNearBy(ivert, layer,
                                                                 theta_center, phi_center, theta_window, phi_window);
    if (Verbosity() >= 1) _t_search_clusters->stop();

#ifdef _DEBUG_
    cout << __LINE__ << ": new_cluster_keys size: " << new_cluster_keys.size() << std::endl;
#endif

    std::vector<PHGenFit::Measurement*> measurements;
    for (TrkrDefs::cluskey cluster_key : new_cluster_keys)
    {
      TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
      if (!cluster)
      {
        LogError("No cluster Found!\n");
        continue;
      }
      
      PHGenFit::Measurement* meas = TrkrClusterToPHGenFitMeasurement(cluster);

      if (meas)
        measurements.push_back(meas);
    }
    std::map<double, shared_ptr<PHGenFit::Track> > incr_chi2s_new_tracks;

#ifdef _DEBUG_
    cout << __LINE__ << ": measurements.size(): " << measurements.size() << endl;
#endif

    if (Verbosity() >= 1) _t_track_propagation->restart();
    track->updateOneMeasurementKalman(measurements, incr_chi2s_new_tracks, extrapolate_base_TP_id, direction, blowup_factor, use_fitted_state);
    use_fitted_state = false;
    blowup_factor = 1.;
    if (Verbosity() >= 1) _t_track_propagation->stop();

#ifdef _DEBUG_
    cout << __LINE__ << ": incr_chi2s_new_tracks.size(): " << incr_chi2s_new_tracks.size() << endl;
#endif

    PHGenFitTrkProp::TrackQuality tq(track_iter->first);

    // Update first track candidate
    if (incr_chi2s_new_tracks.size() > 0)
    {
      auto iter = incr_chi2s_new_tracks.begin();

      if (iter->first < _max_incr_chi2s[layer] and iter->first > 0)
      {
#ifdef _DEBUG_
        cout
            << __LINE__
            << ": iPHGenFitTrack: " << iPHGenFitTrack << endl
            << ": First accepted IncrChi2: " << iter->first << endl
            << "; before update: " << track->get_cluster_keys().back()
            << endl;
#endif
        //				_PHGenFitTracks[iPHGenFitTrack] = std::shared_ptr
        //						< PHGenFit::Track > (iter->second);
        //				track = _PHGenFitTracks[iPHGenFitTrack];

        //				track_iter->first += iter->first;

        track_iter->first.nhits = tq.nhits + 1;
        track_iter->first.chi2 = tq.chi2 + iter->first;
        track_iter->first.ndf = tq.ndf + 2;
        track_iter->first.ntpc = tq.ntpc + (is_tpc_layer(layer) ? 1 : 0);
        track_iter->first.nintt = tq.nintt + (is_intt_layer(layer) ? 1 : 0);
        track_iter->first.nmaps = tq.nmaps + (is_maps_layer(layer) ? 1 : 0);

        track_iter->second = std::shared_ptr<PHGenFit::Track>(iter->second);

        consecutive_missing_layer = 0;
        layer_updated = true;
        extrapolate_base_TP_id = -1;
#ifdef _DEBUG_
        cout
            << __LINE__
            << ": after update: " << track->get_cluster_keys().back()
            << endl;

        fout_chi2
            << _event << "\t"
            << iPHGenFitTrack << "\t"
            << layer << "\t "
            << iter->first
            << endl;
#endif
      }
    }

    // Update other candidates
    if (incr_chi2s_new_tracks.size() > 1 and
        (layer >= _nlayers_maps and layer < _nlayers_maps + _nlayers_intt))
    {
      for (auto iter = (++incr_chi2s_new_tracks.begin());
           iter != incr_chi2s_new_tracks.end(); ++iter)
      {
        if (!(iter->first < _max_splitting_chi2 and iter->first > 0))
          break;

#ifdef _DEBUG_
        std::cout << __LINE__ << ": "
                  << "Track Spliting with "
                  << "IncrChi2: " << iter->first << std::endl;
#endif

        //				_PHGenFitTracks.insert(
        //						MapPHGenFitTrack::value_type(track_iter->first + iter->first,
        //								std::shared_ptr < PHGenFit::Track> (iter->second)));

	// this is a new track, have to associate it with the vertex
	iter->second->set_vertex_id(ivert);
        _PHGenFitTracks.push_back(
            MapPHGenFitTrack::value_type(
                PHGenFitTrkProp::TrackQuality(
                    tq.nhits + 1,
                    tq.chi2 + iter->first,
                    tq.ndf + 2,
                    tq.ntpc + (is_tpc_layer(layer) ? 1 : 0),
                    tq.nintt + (is_intt_layer(layer) ? 1 : 0),
                    tq.nmaps + (is_maps_layer(layer) ? 1 : 0)),
                    iter->second));
      }

#ifdef _DEBUG_
      std::cout << __LINE__ << ": "
                << "_PHGenFitTracksSize: " << _PHGenFitTracks.size() << std::endl;
      std::cout << __LINE__ << ": " << track_iter->second->get_cluster_keys().back() << std::endl;
#endif
    }

#ifdef _DEBUG_
    		cout<<__LINE__<<": updateOneMeasurementKalman:"<<endl;
    		std::cout<<"iPHGenFitTrack: "<<iPHGenFitTrack
    				<<", layer: "<<layer
    				<<", #meas: "<<measurements.size()
    				<<", #tracks: "<<incr_chi2s_new_tracks.size()
		  //<<", #totoal tracks: "<<_trackID_PHGenFitTrack.size()
    				<<std::endl;

    for (auto iter =
             incr_chi2s_new_tracks.begin();
         iter != incr_chi2s_new_tracks.end(); iter++)
    {
      std::cout << __LINE__ << ": IncrChi2: " << iter->first << std::endl;
    }
#endif
    // if (_analyzing_mode)
    // {
    //   int ncand = 0;
    //   for (auto iter =
    //            incr_chi2s_new_tracks.begin();
    //        iter != incr_chi2s_new_tracks.end(); iter++)
    //   {
    //     if (iter->first < _max_incr_chi2s[layer] and iter->first > 0) ncand++;
    //   }
    //   /*
    // 		    float this_pt = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Pt();
    // 		    float this_phi = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Phi();
    // 		    float this_eta = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Eta();
    // 		  */
    //   //"spt:seta:sphi:pt:eta:phi:layer:ncand:nmeas"
    //   //_analyzing_ntuple->Fill(init_pt,init_eta,init_phi,this_pt,this_eta,this_phi,layer,ncand,measurements.size());
    // }
    if (!layer_updated)
      ++consecutive_missing_layer;
  }  // layer loop

#ifdef _DEBUG_
  cout
      << __LINE__
      << ": cluster keys size:  " << track->get_cluster_keys().size()
      << endl;
#endif

  //! Track succesfully propagated and return 0
  return 0;
}

PHGenFit::Measurement* PHGenFitTrkProp::TrkrClusterToPHGenFitMeasurement(
    const TrkrCluster* cluster)
{

  if (!cluster) return nullptr;
  
  TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));
  TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);
  
  // get the trkrid
  TrkrDefs::cluskey cluster_id = cluster->getClusKey();
  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_id);
  int layer = TrkrDefs::getLayer(cluster_id);

  if(trkrid == TrkrDefs::mvtxId)
    {
      int stave_index = MvtxDefs::getStaveId(cluster_id);
      int chip_index = MvtxDefs::getChipId(cluster_id);
      
      double ladder_location[3] = {0.0, 0.0, 0.0};
      CylinderGeom_Mvtx* geom =
	dynamic_cast<CylinderGeom_Mvtx*>(_geom_container_maps->GetLayerGeom(layer));
      // returns the center of the sensor in world coordinates - used to get the ladder phi location
      geom->find_sensor_center(stave_index, 0,
			       0, chip_index, ladder_location);
      //n.Print();
      n.SetXYZ(ladder_location[0], ladder_location[1], 0);
      n.RotateZ(geom->get_stave_phi_tilt());
    }
  else if(trkrid == TrkrDefs::inttId)
    {
      CylinderGeomIntt* geom =
	dynamic_cast<CylinderGeomIntt*>(_geom_container_intt->GetLayerGeom(layer));
      double hit_location[3] = {0.0, 0.0, 0.0};
      geom->find_segment_center(InttDefs::getLadderZId(cluster_id),
				InttDefs::getLadderPhiId(cluster_id), hit_location);
      
      n.SetXYZ(hit_location[0], hit_location[1], 0);
      n.RotateZ(geom->get_strip_phi_tilt());
    }
  
  //double radius = sqrt(cluster->getPosition(0)*cluster->getPosition(0)  + cluster->getPosition(1)*cluster->getPosition(1));
  
  PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
								cluster->getRPhiError(), cluster->getZError());

  meas->set_cluster_key(cluster->getClusKey());

#ifdef _DEBUG_
  int layer_out = TrkrDefs::getLayer(cluster->getClusKey());
  cout
      << __LINE__
      << ": ID: " << cluster->getClusKey()
      << ": layer: " << layer_out
      << ": pos: {" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << "}"
      << ": n: {" << n.X() << ", " << n.Y() << ", " << n.Z() << "}"
      << ": r*phi_error: " << cluster->getRPhiError()
      << ": z error: " << cluster->getZError()
      << endl;
#endif

  return meas;
}

int PHGenFitTrkProp::BuildLayerZPhiHitMap(unsigned int ivert)
{
  // make this map for each collision vertex
  std::multimap<unsigned int, TrkrDefs::cluskey> this_layer_thetaID_phiID_clusterID;

  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;
      
      // This z-phi map relative to the collision vertex includes all clusters, 
      // we do not know whch clusters belong to which vertices yet
      // make a map for every vertex?

      unsigned int layer = TrkrDefs::getLayer(cluskey);
      float x = cluster->getPosition(0) - _vertex[ivert][0];
      float y = cluster->getPosition(1) - _vertex[ivert][1];
      float z = cluster->getPosition(2) - _vertex[ivert][2];
      
      float phi = atan2(y, x);
      float r = sqrt(x * x + y * y);
      float theta = atan2(r, z);

#ifdef _DEBUG_
//		//float rphi = r*phi;
//		std::cout
//		<< __LINE__
//		<<": ID: " << cluster->get_id()
//		<<": layer: "<<cluster->get_layer()
//		<<", r: "<<r
//		<<", z: "<<z
//		<<", phi: "<<phi
//		<<", theta: "<<theta
//		<<endl;
#endif

    unsigned int idx = encode_cluster_index(layer, theta, phi);

// #ifdef _DEBUG_
// 			cout
// 			<<__LINE__<<": "
// 			<<"{ "
// 			<<layer <<", "
// 			<<z <<", "
// 			<<phi << "} =>"
// 			<<idx << ": size: "
// 			<<_layer_thetaID_phiID_clusterID.count(idx)
// 			<<endl;
// #endif

    this_layer_thetaID_phiID_clusterID.insert(std::make_pair(idx, cluster->getClusKey()));
    //cout << "BuildLayerPhiZHitmap for ivert "
    //	 << ivert << " : Inserted pair with  x" << x << " y " << y << " z " << z << " idx " << idx << " cluskey " << cluster->getClusKey() << endl;
  }

  _layer_thetaID_phiID_clusterID.push_back(this_layer_thetaID_phiID_clusterID);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<TrkrDefs::cluskey> PHGenFitTrkProp::SearchHitsNearBy(const unsigned int ivert, const unsigned int layer,
                                                            const float theta_center, const float phi_center, const float theta_window,
                                                            const float phi_window)
{
  if(Verbosity() > 10) cout << "SearchHitsNearBy for ivert = " << ivert << endl;
  
  std::vector<TrkrDefs::cluskey> cluster_keys;
  
  const unsigned int max_phi_bin = 16383;  //2^14 - 1
  const unsigned int max_z_bin = 2047;     // 2^11 - 1
  
  float lower_phi = phi_center - phi_window + _half_max_phi;
  float upper_phi = phi_center + phi_window + _half_max_phi;
  float lower_z = theta_center - theta_window + _half_max_theta;
  float upper_z = theta_center + theta_window + _half_max_theta;

  unsigned int lower_phi_bin = (unsigned int) ((lower_phi) / _layer_thetaID_phiID_clusterID_phiSize);
  unsigned int upper_phi_bin = (unsigned int) ((upper_phi) / _layer_thetaID_phiID_clusterID_phiSize);
  unsigned int lower_z_bin = (unsigned int) ((lower_z) / _layer_thetaID_phiID_clusterID_zSize);
  unsigned int upper_z_bin = (unsigned int) ((upper_z) / _layer_thetaID_phiID_clusterID_zSize);

  if (lower_phi < 0) lower_phi_bin = 0;
  if (upper_phi_bin > max_phi_bin) upper_phi_bin = max_phi_bin;

  if (lower_z < 0) lower_z_bin = 0;
  if (upper_z_bin > max_z_bin) upper_z_bin = max_z_bin;

  for (unsigned int iz = lower_z_bin; iz <= upper_z_bin; ++iz)
  {
    for (unsigned int irphi = lower_phi_bin; irphi <= upper_phi_bin;
         ++irphi)
    {
      if (Verbosity() >= 2) _t_search_clusters_encoding->restart();
      unsigned int idx = encode_cluster_index(layer, iz, irphi);
      if (Verbosity() >= 2) _t_search_clusters_encoding->stop();
      
      if(Verbosity() > 10)
	{
	  if(_layer_thetaID_phiID_clusterID[ivert].count(idx) > 0)
	    {
	      cout
		<<__LINE__<<": "
		<<"{ "
		<<layer <<", "
		<<iz <<", "
		<<irphi << "} =>"
		<<idx << ": size: "
		<<_layer_thetaID_phiID_clusterID[ivert].count(idx)
		<<endl;
	    }
	}
      
      if (Verbosity() >= 2) _t_search_clusters_map_iter->restart();
      for (auto iter = _layer_thetaID_phiID_clusterID[ivert].lower_bound(idx);
      iter != _layer_thetaID_phiID_clusterID[ivert].upper_bound(idx);
      ++iter)
	{
      if(Verbosity() > 10) cout << "      adding cluster with key " << iter->second << endl; 
      cluster_keys.push_back(iter->second);
    }
      if (Verbosity() >= 2) _t_search_clusters_map_iter->stop();
    }
    }
      
      if(Verbosity() > 10)
	{
      //	std::cout
      //			<<__LINE__<<": "
      //			<<"layer: "<<layer
      //			<<": "<<phi_center <<" +- "<<phi_window
      //			<<": "<<theta_center <<" +- "<<theta_window
      //			<<endl;
      std::cout
	<< __LINE__ << ": "
	<< "layer: " << layer
	<< ", rphi: {" << lower_phi_bin
	<< ", " << upper_phi_bin
	<< "}, z: {" << lower_z_bin
	<< ", " << upper_z_bin
	<< "}, found #clusters: " << cluster_keys.size()
	<< endl;
    }

  return cluster_keys;
}

unsigned int PHGenFitTrkProp::encode_cluster_index(const unsigned int layer,
                                                   const float theta, const float phi)
{
  unsigned int idx = UINT_MAX;

  if (layer >= 128)
  {
    LogError("layer >= 128\n");
    return idx;
  }

  if ((theta + _half_max_theta) < 0)
  {
    LogError("(theta + _half_max_theta) < 0 \n");
    return idx;
  }
  unsigned int itheta = (theta + _half_max_theta) / _layer_thetaID_phiID_clusterID_zSize;

  if ((phi + _half_max_phi) < 0)
  {
    LogError("(rphi + _half_max_rphi) < 0 \n");
    return idx;
  }
  unsigned int irphi = (phi + _half_max_phi) / _layer_thetaID_phiID_clusterID_phiSize;

#ifdef _DEBUG_
//	std::cout<<__LINE__<<": "
//			<<": layer: "<<layer
//			<<", irphi: "<<irphi
//			<<", itheta: "<<itheta
//			<<endl;
#endif

  return encode_cluster_index(layer, itheta, irphi);
}

unsigned int PHGenFitTrkProp::encode_cluster_index(const unsigned int layer,
                                                   const unsigned int iz, const unsigned int irphi)
{
  if (layer >= 128)
  {
    LogError("layer >= 128: ") << layer << endl;
    ;
    return UINT_MAX;
  }

  if (iz >= 2048)
  {
    LogError("iz >= 2048: ") << iz << endl;
    return UINT_MAX;
  }

  if (irphi >= 16384)
  {
    LogError("irphi >= 16384: ") << irphi << endl;
    return UINT_MAX;
  }

  unsigned int index = 0;

  index |= (layer << 25);
  index |= (iz << 14);
  index |= (irphi);

  return index;
}

int PHGenFitTrkProp::InitializeGeometry(PHCompositeNode* topNode)
{
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  map<float, int> radius_layer_map;

  _radii_all.assign(_nlayers_all, 0.0);
  _layer_ilayer_map_all.clear();
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange layerrange =
        cellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        mapsladdergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  for (map<float, int>::iterator iter = radius_layer_map.begin();
       iter != radius_layer_map.end(); ++iter)
  {
    _layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));
  }

  // now we extract the information from the cellgeos first
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange begin_end =
        cellgeos->get_begin_end();
    PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++)
    {
      PHG4CylinderCellGeom* geo = miter->second;

      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        laddergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++)
    {
      PHG4CylinderGeom* geo = miter->second;

      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        mapsladdergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++)
    {
      PHG4CylinderGeom* geo = miter->second;

      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

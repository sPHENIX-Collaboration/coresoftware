/*!
 *  \file PHGenFitTrkProp.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#include "PHGenFitTrkProp.h"
#include "AssocInfoContainer.h"

// trackbase_historic includes
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxHit_v1.h>
#include <trackbase_historic/SvtxHitMap.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4Cell.h>
//
#include <g4mvtx/PHG4CylinderGeom_MVTX.h>
#include <g4intt/PHG4CylinderGeomINTT.h>

#include <g4bbc/BbcVertexMap.h>
#include <g4bbc/BbcVertex.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/PHRandomSeed.h>
#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>
 //FIXME remove includes below after having real vertxing
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

// sGeant4 includes
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4FieldManager.hh>

// Helix Hough includes
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/HelixResolution.h>
#include <HelixHough/HelixRange.h>
#include <HelixHough/HelixHough.h>
#include <HelixHough/VertexFinder.h>

// GenFit
#include <GenFit/FieldManager.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>

// gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// standard includes
#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <tuple>
#include <algorithm>
#include <memory>

//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>

#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp

#define _DEBUG_

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
		unsigned int nlayers_tpc
		)
    : PHTrackPropagating(name),
	  _t_seeds_cleanup(nullptr),
	  _t_translate_to_PHGenFitTrack(nullptr),
	  _t_translate1(nullptr),
	  _t_translate2(nullptr),
	  _t_translate3(nullptr),
	  _t_kalman_pat_rec(nullptr),
	  _t_search_clusters(nullptr),
	  _t_search_clusters_encoding(nullptr),
	  _t_search_clusters_map_iter(nullptr),
	  _t_track_propagation(nullptr),
	  _t_full_fitting(nullptr),
	  _t_output_io(nullptr),

      _vertex(),
      _bbc_vertexes(NULL),

			_svtxhitsmap(nullptr),
      _hit_used_map(NULL),
      _cells_svtx(nullptr),
      _cells_intt(nullptr),
      _cells_maps(nullptr),
      _geom_container_intt(nullptr),
      _geom_container_maps(nullptr),
      _analyzing_mode(false),
      _analyzing_file(NULL),
      _analyzing_ntuple(NULL),
      _max_merging_dphi(0.1),
      _max_merging_deta(0.1),
      _max_merging_dr(0.1),
      _max_merging_dz(0.1),
      _max_share_hits(3),
      _fitter(NULL),
      _track_fitting_alg_name("KalmanFitter"),
      _primary_pid_guess(211),
      _cut_min_pT(0.2),
      _do_evt_display(false),
      _nlayers_maps(nlayers_maps),
      _nlayers_intt(nlayers_intt),
      _nlayers_tpc(nlayers_tpc),
      _nlayers_all(_nlayers_maps+_nlayers_intt+_nlayers_tpc),
      _layer_ilayer_map_all(),
      _radii_all(),
      
      _max_search_win_phi_tpc(    0.0040),
      _min_search_win_phi_tpc(    0.0000),
      _max_search_win_theta_tpc(  0.0040),
      _min_search_win_theta_tpc(  0.0000),      
      
      _max_search_win_phi_maps(   0.0050),
      _min_search_win_phi_maps(   0.0000),
      _max_search_win_theta_maps( 0.0400),
      _min_search_win_theta_maps( 0.0000),
      
      _search_win_phi(20),
      _search_win_theta(20),
      _layer_thetaID_phiID_cluserID(),
      //_half_max_theta(160),
      _half_max_theta(3.1416/2.),
      //_half_max_phi(252), //80cm * Pi
      _half_max_phi(3.1416),
      //_layer_thetaID_phiID_cluserID_phiSize(0.1200),
      _layer_thetaID_phiID_cluserID_phiSize(0.1200/30), //rad
      _layer_thetaID_phiID_cluserID_zSize(0.1700/30),
      _PHGenFitTracks(),
      _init_direction(-1),
      _blowup_factor(1.),
      _max_consecutive_missing_layer(20),
      _max_incr_chi2(20.),
      _max_splitting_chi2(20.),
      _min_good_track_hits(30)
{
  _event = 0;

	_max_search_win_phi_intt[0] = 0.20;
	_max_search_win_phi_intt[1] = 0.20;
	_max_search_win_phi_intt[2] = 0.0050;
	_max_search_win_phi_intt[3] = 0.0050;
	_max_search_win_phi_intt[4] = 0.0050;
	_max_search_win_phi_intt[5] = 0.0050;
	_max_search_win_phi_intt[6] = 0.0050;
	_max_search_win_phi_intt[7] = 0.0050;

	_min_search_win_phi_intt[0] =   0.2000;
	_min_search_win_phi_intt[1] =   0.2000;
	_min_search_win_phi_intt[2] =   0.0;
	_min_search_win_phi_intt[3] =   0.0;
	_min_search_win_phi_intt[4] =   0.0;
	_min_search_win_phi_intt[5] =   0.0;
	_min_search_win_phi_intt[6] =   0.0;
	_min_search_win_phi_intt[7] =   0.0;

	_max_search_win_theta_intt[0] = 0.010;
	_max_search_win_theta_intt[1] = 0.010;
	_max_search_win_theta_intt[2] =  0.2000;
	_max_search_win_theta_intt[3] =  0.2000;
	_max_search_win_theta_intt[4] =  0.2000;
	_max_search_win_theta_intt[5] =  0.2000;
	_max_search_win_theta_intt[6] =  0.2000;
	_max_search_win_theta_intt[7] =  0.2000;

	_min_search_win_theta_intt[0] = 0.000;
	_min_search_win_theta_intt[1] = 0.000;
	_min_search_win_theta_intt[2] = 0.200;
	_min_search_win_theta_intt[3] = 0.200;
	_min_search_win_theta_intt[4] = 0.200;
	_min_search_win_theta_intt[5] = 0.200;
	_min_search_win_theta_intt[6] = 0.200;
	_min_search_win_theta_intt[7] = 0.200;

	_vertex_error.clear();
	_vertex_error.assign(3, 0.0100);
}

int PHGenFitTrkProp::Setup(PHCompositeNode* topNode) {

	// Start new interface ----
	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = PHTrackPropagating::Setup(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = CreateNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;
	// End new interface ----

	ret = InitializePHGenFit(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK)
		return ret;

	if(_analyzing_mode){
	  cout << "Ana Mode, creating ntuples! " << endl;
	  _analyzing_file = new TFile("./PHGenFitTrkProp.root","RECREATE");
	  //	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","spt:seta:sphi:pt:eta:phi:layer:ncand:nmeas");
	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","pt:kappa:d:phi:dzdl:z0:nhit:ml:rec:dt");
	  cout << "Done" << endl;
	  
	}
	
//	ret = InitializeGeometry(topNode);
//	if(ret != Fun4AllReturnCodes::EVENT_OK)
//	  return ret;

	/*!
	 * Initilize parameters
	 */
	for(int layer = 0; layer < _nlayers_all; ++layer) {
		_search_wins_phi.insert(std::make_pair(layer, _search_win_phi));
		_search_wins_theta.insert(std::make_pair(layer, _search_win_theta));
		_max_incr_chi2s.insert(std::make_pair(layer, _max_incr_chi2));
	}

#ifdef _DEBUG_
	for(int layer = 0; layer < _nlayers_all; ++layer) {
		cout
		<<__LINE__
		<<": layer: "<< layer
		<<": search_wins_rphi: " << _search_wins_phi[layer]
		<<": search_wins_z: " << _search_wins_theta[layer]
		<<": max_incr_chi2: " << _max_incr_chi2s[layer]
		<<endl;
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

int PHGenFitTrkProp::Process() {

  if (Verbosity() > 0){
	  cout << "PHGenFitTrkProp::process_event -- entered" << endl;
	  cout << "nMapsLayers = " << _nlayers_maps << endl;
	  cout << "nInttLayers = " << _nlayers_intt << endl;
	  cout << "nTPCLayers = " << _nlayers_tpc << endl;
  }
	// start fresh
	int ret;

	// TODO vertex using strategy
	_vertex.clear();
	_vertex.assign(3, 0.0);

	if(_vertex_map){
		SvtxVertex* vertex = _vertex_map->get(0);
		TVector3 v(vertex->get_x(),vertex->get_y(),vertex->get_z());
		_vertex[0] = vertex->get_x();
		_vertex[1] = vertex->get_y();
		_vertex[2] = vertex->get_z();
		for(int i=0;i<3;++i)
			_vertex_error[i] = sqrt(vertex->get_error(i,i));
	}

	{
	  //-----------------------------------
	  // Kalman track propagating
	  //-----------------------------------
	  if(Verbosity() >= 1) _t_kalman_pat_rec->restart();
		ret = KalmanTrkProp();
		if (ret != Fun4AllReturnCodes::EVENT_OK)
			return ret;
	  if(Verbosity() >= 1) _t_kalman_pat_rec->stop();
	}

  if(Verbosity() > 1) print_timers();

	++_event;
	
	return Fun4AllReturnCodes::EVENT_OK;
}

void PHGenFitTrkProp::print_timers() {
  
  std::cout << "=============== Timers: ===============" << std::endl;
  std::cout << "\t - Seeds Cleanup:          "<<_t_seeds_cleanup->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "CPUSCALE Pattern recognition time:    "<<_t_kalman_pat_rec->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t - Track Translation time: "<<_t_translate_to_PHGenFitTrack->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t -    - Translation1 time: "<<_t_translate1->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t -    - Translation2 time: "<<_t_translate2->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t -    - Translation3 time: "<<_t_translate3->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t - Cluster searching time: "<<_t_search_clusters->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t\t - Encoding time:        "<<_t_search_clusters_encoding->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t\t - Map iteration:        "<<_t_search_clusters_map_iter->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "\t - Kalman updater time:    "<<_t_track_propagation->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "Full fitting time:           "<<_t_full_fitting->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "Output IO time:              "<<_t_output_io->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "=======================================" << std::endl;

}

int PHGenFitTrkProp::End() {

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

	if(_analyzing_mode){
	  cout << " cleaning up " << endl;
	  _analyzing_file->cd();
	  _analyzing_ntuple->Write();
	  _analyzing_file->Close();
	  //	  delete _analyzing_ntuple;
	  // delete _analyzing_file;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::CreateNodes(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::InitializePHGenFit(PHCompositeNode* topNode) {

  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);


  PHField * field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager, field, _track_fitting_alg_name,
					"RKTrackRep", _do_evt_display);

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

#ifdef _DEBUG_
	_fitter->set_verbosity(10);
#endif

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::GetNodes(PHCompositeNode* topNode) {

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	// used in fast vertexing from BBC
	_bbc_vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");

  // get node containing the digitized hits
  _svtxhitsmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_svtxhitsmap) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _cells_svtx = findNode::getClass<PHG4CellContainer>(topNode,
						      "G4CELL_TPC");

  _cells_intt = findNode::getClass<PHG4CellContainer>(
						      topNode, "G4CELL_INTT");

  _cells_maps = findNode::getClass<PHG4CellContainer>(
						      topNode, "G4CELL_MVTX");

  if (!_cells_svtx and !_cells_intt and !_cells_maps) {
    if (Verbosity() >= 0) {
      LogError("No PHG4CellContainer found!");}
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _geom_container_intt = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  _geom_container_maps = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  if (!_cells_svtx && !_cells_maps && !_cells_intt) {
    cout << PHWHERE << "ERROR: Can't find any cell node!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::check_track_exists(MapPHGenFitTrack::iterator iter){
	
  
  //Loop over hitIDs on current track and check if they have been used
  unsigned int n_clu = iter->second->get_cluster_IDs().size();

  unsigned int  n_clu_used = 0;
  const std::vector<unsigned int>& clusterIDs = iter->second->get_cluster_IDs();
  for(unsigned int iCluId = 0; iCluId < clusterIDs.size(); ++iCluId){
    unsigned int cluster_ID = clusterIDs[iCluId];
    //if(_hit_used_map[cluster_ID]>0) n_clu_used++;
    if(_assoc_container->GetTracksFromCluster(cluster_ID).size()>0) n_clu_used++;
  }
  int code = 0;
  if(((float)n_clu_used/n_clu)>0.3){
    if(Verbosity()>=1)
      cout << "Found duplicate track. n_clu: " << n_clu << " c_clu_used: " << n_clu_used << endl;
    /*
    for(unsigned int iCluId = 0; iCluId < clusterIDs.size(); ++iCluId){
      unsigned int cluster_ID = clusterIDs[iCluId];
      cout << "#Clu_g = " << iCluId 
	   << " layer: " << _cluster_map->get(cluster_ID)->get_layer()
	   << " r: " << TMath::Sqrt(_cluster_map->get(cluster_ID)->get_x()*_cluster_map->get(cluster_ID)->get_x() +_cluster_map->get(cluster_ID)->get_y()*_cluster_map->get(cluster_ID)->get_y() )
	   << endl;
    }
    */
    return code;
  }
  code = 1;
  return code;
}

int PHGenFitTrkProp::KalmanTrkProp() {

#ifdef _DEBUG_
	std::cout << "=========================" << std::endl;
	std::cout << "PHGenFitTrkProp::KalmanTrkProp: Start: Event: "<< _event << std::endl;
	std::cout << "Total Raw Tracks: " <<_track_map->size() << std::endl;
	std::cout << "=========================" << std::endl;
#endif

	/*!
	 *   sort clusters
	 */
	BuildLayerZPhiHitMap();

	vector<genfit::Track*> evt_disp_copy;

	_PHGenFitTracks.clear();


	for (auto trackiter = _track_map->begin();
			trackiter != _track_map->end(); ++trackiter) {

		SvtxTrack *tracklet = trackiter->second;

#ifdef _DEBUG_
		std::cout
		<< __LINE__
		<< ": Processing itrack: " << trackiter->first
		<< ": Total tracks: " << _track_map->size()
		<<endl;
#endif

		/*!
		 * Translate SimpleTrack3D To PHGenFitTracks
		 */
		if(Verbosity() > 1) _t_translate_to_PHGenFitTrack->restart();
		SvtxTrackToPHGenFitTracks(tracklet);
		if(Verbosity() > 1) _t_translate_to_PHGenFitTrack->stop();

		/*!
		 * Handle track propagation, termination, output and evt disp.
		 */
		bool is_splitting_track = false;
#ifdef _DEBUG_
		int i = 0;
#endif
		for (auto iter = _PHGenFitTracks.begin();
				iter != _PHGenFitTracks.end(); ++iter) {
#ifdef _DEBUG_
			cout
			<< __LINE__
			<< ": propergating: " << i <<"/" << trackiter->first
			<< endl;

			cout << _cluster_map << endl;
			if(_cluster_map) _cluster_map->identify();
			else cout << __LINE__ << "_cluster_map not init!" << endl;
#endif

			std::vector<unsigned int> clusterIDs = iter->second->get_cluster_IDs();

			unsigned int init_layer = UINT_MAX;

			if(!is_splitting_track) {
				if(_init_direction == 1) {
					init_layer = _cluster_map->get(clusterIDs.front())->get_layer();
					TrackPropPatRec(iter, init_layer, _nlayers_all, true);
					TrackPropPatRec(iter, init_layer, 0, false);
				} else {
					init_layer = _cluster_map->get(clusterIDs.back())->get_layer();
					TrackPropPatRec(iter, init_layer, 0, true);
					TrackPropPatRec(iter, init_layer, _nlayers_all, false);
				}
				is_splitting_track = true;
			} else {
				if(_init_direction == 1) {
					init_layer = _cluster_map->get(clusterIDs.front())->get_layer();
					TrackPropPatRec(iter, init_layer, _nlayers_all, false);
				} else {
					init_layer = _cluster_map->get(clusterIDs.back())->get_layer();
					TrackPropPatRec(iter, init_layer, 0, false);
				}
			}

#ifdef _DEBUG_
			cout
			<< __LINE__
			<< ": tracki: " << i
			<< ": clusterIDs size:  " << iter->second->get_cluster_IDs().size()
			<< ": quality: " << iter->first
			<< endl;
			++i;
#endif

			//_trackID_PHGenFitTrack.erase(iter);
		}// loop _PHGenFitTracks

		if(_PHGenFitTracks.size()==0) continue;

#ifdef _DEBUG_
		i = 0;
		for (auto iter = _PHGenFitTracks.begin();
				iter != _PHGenFitTracks.end(); ++iter) {
			cout
			<< __LINE__
			<< ": track: " << i++
			<< ": clusterIDs size:  " << iter->second->get_cluster_IDs().size()
			<< ": quality: " << iter->first
			<< endl;
		}
#endif

		//std::sort(_PHGenFitTracks.begin(), _PHGenFitTracks.end());
		_PHGenFitTracks.sort();

#ifdef _DEBUG_
		for (auto iter = _PHGenFitTracks.begin();
				iter != _PHGenFitTracks.end(); ++iter) {
			cout
			<< __LINE__
			<< ": clusterIDs size:  " << iter->second->get_cluster_IDs().size()
			<< ": quality: " << iter->first
			<< endl;
		}
#endif

		auto iter = _PHGenFitTracks.begin();
		
		int track_exists = check_track_exists(iter);
		if (iter->second->get_cluster_IDs().size() >= _min_good_track_hits && track_exists) {
			OutputPHGenFitTrack(iter);
#ifdef _DEBUG_
			cout << __LINE__ << endl;
#endif
			if (_do_evt_display) {
				evt_disp_copy.push_back(
						new genfit::Track(*iter->second->getGenFitTrack()));
			}
		}

		_PHGenFitTracks.clear();
	}

#ifdef _DEBUG_
	std::cout << "=========================" << std::endl;
	std::cout << "PHGenFitTrkProp::KalmanTrkProp: End: Event: "<< _event << std::endl;
	std::cout << "Total Final Tracks: " << _track_map->size() << std::endl;
	std::cout << "=========================" << std::endl;
#endif

	if (_do_evt_display) {
		_fitter->getEventDisplay()->addEvent(evt_disp_copy);
	} else {
		evt_disp_copy.clear();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::OutputPHGenFitTrack(MapPHGenFitTrack::iterator iter) {

//	for (std::map<int, std::shared_ptr<PHGenFit::Track>>::iterator iter =
//			_trackID_PHGenFitTrack.begin();
//			iter != _trackID_PHGenFitTrack.end(); iter++) {

#ifdef _DEBUG_
		std::cout << "=========================" << std::endl;
		//std::cout << __LINE__ << ": iPHGenFitTrack: " << iter->first << std::endl;
		std::cout << __LINE__ << ": _track_map->size(): " << _track_map->size() << std::endl;
		std::cout << "Contains: " << iter->second->get_cluster_IDs().size() << " clusters." <<std::endl;
		std::cout << "=========================" << std::endl;
#endif

		SvtxTrack_v1 track;
		//track.set_id(iter->first);
		track.set_id(_track_map->size());

#ifdef _DO_FULL_FITTING_
		if(Verbosity() >= 1) _t_full_fitting->restart();
		if (_fitter->processTrack(iter->second.get(), false) != 0) {
			if (Verbosity() >= 1)
				LogWarning("Track fitting failed\n");
			//delete track;
			return -1;
		}
		if(Verbosity() >= 1) _t_full_fitting->stop();

		if(Verbosity() >= 1) _t_output_io->restart();
//		iter->second->getGenFitTrack()->Print();

		track.set_chisq(iter->second->get_chi2());
		track.set_ndf(iter->second->get_ndf());

		//FIXME use fitted vertex
		TVector3 vertex_position(0, 0, 0);
		std::unique_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca = NULL;
		try {
			gf_state_vertex_ca = std::unique_ptr < genfit::MeasuredStateOnPlane
					> (iter->second->extrapolateToPoint(vertex_position));
		} catch (...) {
			if (Verbosity() >= 2)
				LogWarning("extrapolateToPoint failed!");
		}
		if (!gf_state_vertex_ca) {
			//delete out_track;
			return -1;
		}

		TVector3 mom = gf_state_vertex_ca->getMom();
		TVector3 pos = gf_state_vertex_ca->getPos();
		TMatrixDSym cov = gf_state_vertex_ca->get6DCov();
#else
		TVectorD state = iter->second->getGenFitTrack()->getStateSeed();
		TVector3 pos(state(0),state(1),state(2));
		TVector3 mom(state(3),state(4),state(5));
#endif
		track.set_px(mom.Px());
		track.set_py(mom.Py());
		track.set_pz(mom.Pz());

		track.set_x(pos.X());
		track.set_y(pos.Y());
		track.set_z(pos.Z());

		for(unsigned int cluster_ID : iter->second->get_cluster_IDs()){
		  track.insert_cluster(cluster_ID);
		}

		//Check track quality
		//		bool is_good_track = true;

		Int_t n_maps = 0;
		Int_t n_intt = 0;
		Int_t n_tpc  = 0;
		
		for (SvtxTrack::ConstClusterIter iter = track.begin_clusters();
		     iter != track.end_clusters();
		     ++iter) {
		  unsigned int cluster_id = *iter;
		  SvtxCluster* cluster = _cluster_map->get(cluster_id);
		  unsigned int layer = cluster->get_layer();
		  if(_nlayers_maps>0&&layer<_nlayers_maps){ 
		    n_maps++ ;
		  }
		  if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){
		    n_intt++;
		  }
		  if(n_intt >8)
		    {
		      cout << PHWHERE << " Can not have more than 8 INTT layers, quit!" << endl;
		      exit(1);
		    }
		  if(_nlayers_tpc>0&&
		     layer>=(_nlayers_maps+_nlayers_intt)&&
		     layer<(_nlayers_maps+_nlayers_intt+_nlayers_tpc)){ 
		    n_tpc++;
		  }
		}
		/*
		  if(n_maps<3&&_nlayers_maps>0) is_good_track = false;
		  if(n_intt<3&&_nlayers_intt>0) is_good_track = false;
		  if(n_tpc<20&&_nlayers_tpc>0) is_good_track = false;
		*/	
		//if(is_good_track||_n_iteration==4)
		//if(is_good_track||_n_iteration>=0)
		//if(_n_iteration>=0)
		{
			for(unsigned int cluster_ID : iter->second->get_cluster_IDs()){
				//_hit_used_map[cluster_ID] = _n_iteration;
				_assoc_container->SetClusterTrackAssoc(cluster_ID, track.get_id());
			}

			_track_map->insert(&track);
		}
		if (Verbosity() > 5) {
			cout << "track " << _track_map->size() << " quality = " << track.get_quality()
					<< endl;
			cout << "px = " << track.get_px() << " py = " << track.get_py()
					<< " pz = " << track.get_pz() << endl;
		}

		if(Verbosity() >= 1) _t_output_io->stop();
//	}



	return 0;
}

int PHGenFitTrkProp::SvtxTrackToPHGenFitTracks(const SvtxTrack* svtxtrack) {

	// clean up working array for each event
	_PHGenFitTracks.clear();

	double time1 = 0;
	double time2 = 0;

	//FIXME used in analysis ntuple
	float kappa = 0;
	float d = 0;
	float phi = 0;
	float dzdl = 0;
	float z0 = 0;
	float nhit = 0;
	float ml = 0;
	float rec = 0;
	float dt = 0;

	if(Verbosity() > 1){
	  time1 = _t_translate1->get_accumulated_time();
	  time2 = _t_translate1->get_accumulated_time();
	  _t_translate1->restart();
	}

	TVector3 seed_pos(
			svtxtrack->get_x(),
			svtxtrack->get_y(),
			svtxtrack->get_z()
			);

	TVector3 seed_mom(
			svtxtrack->get_px(),
			svtxtrack->get_py(),
			svtxtrack->get_pz()
			);

	const float blowup_factor = 1.;

	TMatrixDSym seed_cov(6);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			seed_cov[i][j] = blowup_factor*svtxtrack->get_error(i, j);
		}
	}

	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
	std::shared_ptr<PHGenFit::Track> track(
			new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov));

	std::multimap<float, unsigned int> m_r_clusterID;
	for (auto hit_iter = svtxtrack->begin_clusters();
			hit_iter!= svtxtrack->end_clusters(); ++hit_iter) {

		SvtxCluster *cluster = _cluster_map->get(*hit_iter);

		float r = sqrt(
				cluster->get_x() * cluster->get_x() +
				cluster->get_y() * cluster->get_y());

		m_r_clusterID.insert(std::pair<float, unsigned int>(r, cluster->get_id()));
	}

	std::vector<PHGenFit::Measurement*> measurements;
	if(_vertex_map){
		TVector3 v(_vertex[0],_vertex[1],_vertex[2]);
		TMatrixDSym cov(3);
		cov.Zero();
		cov(0, 0) = _vertex_error[0]*_vertex_error[0];
		cov(1, 1) = _vertex_error[1]*_vertex_error[1];
		cov(2, 2) = _vertex_error[2]*_vertex_error[2];
		PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(v, cov);
		//FIXME re-use the first cluster id
		unsigned int id = m_r_clusterID.begin()->second;
		meas->set_cluster_ID(id);
		measurements.push_back(meas);
	}

	for (auto iter = m_r_clusterID.begin();
			iter != m_r_clusterID.end();
			++iter) {

		unsigned int cluster_ID = iter->second;

		SvtxCluster* cluster = _cluster_map->get(cluster_ID);
		ml+=cluster->get_layer();
		if (!cluster) {
			LogError("No cluster Found!\n");continue;
		}

		PHGenFit::Measurement *meas = SvtxClusterToPHGenFitMeasurement(cluster);

		if(meas)
			measurements.push_back(meas);
	}
	track->addMeasurements(measurements);

	if(Verbosity() > 1) _t_translate2->stop();
	if(Verbosity() > 1) _t_translate3->restart();

	if (_fitter->processTrack(track.get(), false) != 0) {
		if (Verbosity() >= 1)
			LogWarning("Seed fitting failed")<<std::endl;
		if(Verbosity() > 1) _t_translate3->stop();
		if(Verbosity() > 1){
		  _t_translate1->stop();
		  time2 = _t_translate1->get_accumulated_time();
		}
		dt = time2 - time1;
		if(_analyzing_mode == true)
		  _analyzing_ntuple->Fill(svtxtrack->get_pt(),kappa,d,phi,dzdl,z0,nhit,ml/nhit,rec,dt);
		return -1;
	}

	int nhits = track->get_cluster_IDs().size();
	float chi2 = track->get_chi2();
	float ndf  = track->get_ndf();

	if(nhits > 0 and chi2 > 0 and ndf > 0) {
		_PHGenFitTracks.push_back(
				MapPHGenFitTrack::value_type(
						PHGenFitTrkProp::TrackQuality(nhits, chi2, ndf, nhits, 0, 0), track)
		);
	}
	if(Verbosity() > 1) _t_translate3->stop();
	if(Verbosity() > 1){
	  _t_translate1->stop();
	  time2 = _t_translate1->get_accumulated_time();
	}
	dt = time2 - time1;
	rec = 1;
	if(_analyzing_mode == true)
	  _analyzing_ntuple->Fill(svtxtrack->get_pt(),kappa,d,phi,dzdl,z0,nhit,rec,dt);

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHGenFitTrkProp::TrackPropPatRec(
		MapPHGenFitTrack::iterator &track_iter,
		unsigned int init_layer, unsigned int end_layer,
		const bool use_fitted_state_once) {

#ifdef _DEBUG_
	cout
	<< __LINE__
	<< " TrackPropPatRec"
	<< " : init_layer: " << init_layer
	<< " : end_layer: " << end_layer
	<< " : use_fitted_state_once: " << use_fitted_state_once
	<< endl;
#endif

	std::shared_ptr<PHGenFit::Track> &track = track_iter->second;

	int direction = end_layer >= init_layer ? 1 : -1;
	assert(direction==1 or direction==-1);

	int first_extrapolate_base_TP_id = -1;

	bool use_fitted_state = use_fitted_state_once;
	float blowup_factor = use_fitted_state? _blowup_factor : 1.;

	/*!
	 * Find the last layer of with TrackPoint (TP)
	 * Asumming measuremnts are sorted by radius
	 * and cluster IDs are syncronized with the TP IDs
	 */
	{
		std::vector<unsigned int> clusterIDs = track->get_cluster_IDs();

		for (unsigned int i = 0; i < clusterIDs.size(); ++i) {
			if (_cluster_map->get(clusterIDs[i])->get_layer() == init_layer) {
				first_extrapolate_base_TP_id = i;
				break;
			}
		}
	}

	if(first_extrapolate_base_TP_id < 0) {
		if(Verbosity() > 0)
			LogError("first_extrapolate_base_TP_id < 0");
		return -1;
	}


	int extrapolate_base_TP_id = first_extrapolate_base_TP_id ;

	unsigned int consecutive_missing_layer = 0;

//	unsigned int layer = init_layer + direction;
//	while (layer>=0 and layer < (unsigned int)_nlayers_all and layer!=end_layer) {

	for(unsigned int layer = init_layer + direction;
			layer != end_layer+direction;
			layer += direction) {
		if(!(layer>=0 and layer < (unsigned int)_nlayers_all)) break;

//		if(layer >= _nlayers_maps and layer < _nlayers_maps+_nlayers_intt) continue;

		/*!
		 * if miss too many layers terminate track propagating
		 */
		if(consecutive_missing_layer > _max_consecutive_missing_layer) {
			if(Verbosity() > 1) {
				LogWarning ("consecutive_missing_layer > ") << _max_consecutive_missing_layer << endl;
			}
			if(track->get_cluster_IDs().size() >= _min_good_track_hits)
				return 0;
			else
				return -1;
		}

		bool layer_updated = false;

		float layer_r = _radii_all[_layer_ilayer_map_all[layer]];

#ifdef _DEBUG_
		const int iPHGenFitTrack = _PHGenFitTracks.size();
		std::cout<<"========================="<<std::endl;
		std::cout<<__LINE__<<": Event: "<< _event <<": _PHGenFitTracks.size(): "<<_PHGenFitTracks.size() <<": layer: "<<layer<<std::endl;
		std::cout<<"========================="<<std::endl;
#endif

#ifdef _DEBUG_
		{
			unsigned int tempIdx = extrapolate_base_TP_id >= 0 ? extrapolate_base_TP_id : extrapolate_base_TP_id + track->get_cluster_IDs().size();
			cout
			<< __LINE__
			<<" tempIdx: " <<tempIdx
			<<endl;
			if(tempIdx>=0 and tempIdx < track->get_cluster_IDs().size()) {
				unsigned int extrapolate_base_cluster_id = track->get_cluster_IDs()[tempIdx];
				SvtxCluster* extrapolate_base_cluster = _cluster_map->get(extrapolate_base_cluster_id);
				cout
				<<__LINE__
				<<": Target layer: { " << layer
				<<", " << layer_r
				<<"} : From layer: { " << extrapolate_base_cluster->get_layer()
				<<", " << sqrt(extrapolate_base_cluster->get_x()*extrapolate_base_cluster->get_x() + extrapolate_base_cluster->get_y()*extrapolate_base_cluster->get_y())
				<<"} : ID: " << extrapolate_base_cluster_id
				<<endl;
			}
		}
#endif

//		bool have_tp_with_fit_info = false;
//		std::vector<unsigned int> clusterIDs = track->get_cluster_IDs();
//		for (unsigned int i = clusterIDs.size() - 1; i >= 0; --i) {
//			std::unique_ptr<genfit::MeasuredStateOnPlane> kfsop = NULL;
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
		try {
			state = std::unique_ptr < genfit::MeasuredStateOnPlane
					> (track->extrapolateToCylinder(layer_r, TVector3(0, 0, 0),
							TVector3(0, 0, 1), extrapolate_base_TP_id, direction));
//		genfit::MeasuredStateOnPlane *state = track->extrapolateToCylinder(
//				layer_r, TVector3(0, 0, 0), TVector3(0, 0, 1), 0);
		} catch (...) {
			if (Verbosity() > 1) {
				LogWarning("Can not extrapolate to Cylinder!")<<std::endl;
			}
			continue;
		}

		if(!state) {
			if (Verbosity() > 1) {
				LogWarning("Can not extrapolate to Cylinder!")<<std::endl;
			}
			continue;
		}

		TVector3 pos = state->getPos();
		pos.SetXYZ(
				pos.X()-_vertex[0],
				pos.Y()-_vertex[1],
				pos.Z()-_vertex[2]
				);

		float phi_center = pos.Phi();
		float theta_center = pos.Theta();

#ifdef _USE_CONSTANT_SEARCH_WIN_

		float phi_window     = 25e-4;
		float theta_window   = 25e-4;

		if(layer >= 3 and layer <=6) {
			phi_window     = 300e-4;
			theta_window   = 0.2;
		}

		if(layer <=2 ){
			phi_window     = 3000e-4;
			theta_window   = 3000e-4;
		}
#else
		TMatrixDSym cov = state->get6DCov();

		float phi_window     = _search_wins_phi[layer] * sqrt(cov[0][0] + cov[1][1] + cov[0][1] + cov[1][0]) / pos.Perp();
		float theta_window   = _search_wins_theta[layer]    * sqrt(cov[2][2]) / pos.Perp();

		if(layer < _nlayers_maps){
			if (phi_window > _max_search_win_phi_maps) phi_window = _max_search_win_phi_maps;
			if (phi_window < _min_search_win_phi_maps) phi_window = _min_search_win_phi_maps;
			if (theta_window   > _max_search_win_theta_maps)   theta_window   = _max_search_win_theta_maps;
			if (theta_window   < _min_search_win_theta_maps)   theta_window   = _min_search_win_theta_maps;
		} else if(layer < _nlayers_maps + _nlayers_intt) {
			if (phi_window > _max_search_win_phi_intt[layer - _nlayers_maps]) phi_window = _max_search_win_phi_intt[layer - _nlayers_maps];
			if (phi_window < _min_search_win_phi_intt[layer - _nlayers_maps]) phi_window = _min_search_win_phi_intt[layer - _nlayers_maps];
			if (theta_window   > _max_search_win_theta_intt[layer - _nlayers_maps])   theta_window   = _max_search_win_theta_intt[layer - _nlayers_maps];
			if (theta_window   < _min_search_win_theta_intt[layer - _nlayers_maps])   theta_window   = _min_search_win_theta_intt[layer - _nlayers_maps];
		} else {
			if (phi_window > _max_search_win_phi_tpc) phi_window = _max_search_win_phi_tpc;
			if (phi_window < _min_search_win_phi_tpc) phi_window = _min_search_win_phi_tpc;
			if (theta_window   > _max_search_win_theta_tpc)   theta_window   = _max_search_win_theta_tpc;
			if (theta_window   < _min_search_win_theta_tpc)   theta_window   = _min_search_win_theta_tpc;
		}

		//FIXME optimize this
//		if(layer == _nlayers_maps + _nlayers_intt -1) {
//			phi_window = 0.02;
//			theta_window = 0.04;
//		}

#endif

#ifdef _DEBUG_
		cout<<__LINE__<<": ";
		printf("layer: %d: r: %f: phi: %f +- %f; theta: %f +- %f\n",
				layer, pos.Perp(),
				phi_center, phi_window,
				theta_center, theta_window
				);

//		cout<<__LINE__<<": ";
//		printf("layer: %d:  phi: %f +- %f\n",
//				layer,
//				pos.Phi(), phi_window
//				);
#endif

		if(Verbosity() >= 1) _t_search_clusters->restart();
		std::vector<unsigned int> new_cluster_IDs = SearchHitsNearBy(layer,
				theta_center, phi_center, theta_window, phi_window);
		if(Verbosity() >= 1) _t_search_clusters->stop();

#ifdef _DEBUG_
		cout<<__LINE__<< ": new_cluster_IDs size: " << new_cluster_IDs.size() << std::endl;
#endif

		std::vector<PHGenFit::Measurement*> measurements;
		for (unsigned int cluster_ID : new_cluster_IDs) {
			//LogDebug("cluster_ID: ")<<cluster_ID<<endl;
			SvtxCluster* cluster = _cluster_map->get(cluster_ID);
			if (!cluster) {
				LogError("No cluster Found!\n");
				continue;
			}

			PHGenFit::Measurement *meas = SvtxClusterToPHGenFitMeasurement(cluster);

			if(meas)
				measurements.push_back(meas);
		}
		//std::map<double, PHGenFit::Track*> incr_chi2s_new_tracks;
		std::map<double, shared_ptr<PHGenFit::Track> > incr_chi2s_new_tracks;

#ifdef _DEBUG_
		cout<<__LINE__<<": measurements.size(): "<<measurements.size()<<endl;
#endif

		if(Verbosity() >= 1) _t_track_propagation->restart();
		track->updateOneMeasurementKalman(measurements, incr_chi2s_new_tracks, extrapolate_base_TP_id, direction, blowup_factor, use_fitted_state);
		use_fitted_state = false;
		blowup_factor = 1.;
		if(Verbosity() >= 1) _t_track_propagation->stop();

#ifdef _DEBUG_
		cout<<__LINE__<<": incr_chi2s_new_tracks.size(): "<<incr_chi2s_new_tracks.size()<<endl;
#endif
		
		PHGenFitTrkProp::TrackQuality tq(track_iter->first);

		// Update first track candidate
		if (incr_chi2s_new_tracks.size() > 0) {
			auto iter = incr_chi2s_new_tracks.begin();

			if (iter->first < _max_incr_chi2s[layer] and iter->first > 0) {

#ifdef _DEBUG_
				cout
				<< __LINE__
				<< ": iPHGenFitTrack: " << iPHGenFitTrack << endl
				<< ": First accepted IncrChi2: " << iter->first << endl
				<< "; before update: " << track->get_cluster_IDs().back()
				<< endl;
#endif
//				_PHGenFitTracks[iPHGenFitTrack] = std::shared_ptr
//						< PHGenFit::Track > (iter->second);
//				track = _PHGenFitTracks[iPHGenFitTrack];

//				track_iter->first += iter->first;

				track_iter->first.nhits = tq.nhits + 1;
				track_iter->first.chi2  = tq.chi2  + iter->first;
				track_iter->first.ndf   = tq.ndf   + 2;
				track_iter->first.ntpc  = tq.ntpc  + ((layer >= _nlayers_maps + _nlayers_intt) ? 1 : 0);
				track_iter->first.nintt = tq.nintt + ((layer >= _nlayers_maps and layer < _nlayers_maps + _nlayers_intt) ? 1 : 0);
				track_iter->first.nmaps = tq.nmaps + ((layer < _nlayers_maps) ? 1 : 0);

				track_iter->second = std::shared_ptr<PHGenFit::Track> (iter->second);

				consecutive_missing_layer = 0;
				layer_updated = true;
				extrapolate_base_TP_id = -1;
#ifdef _DEBUG_
				cout
				<< __LINE__
				<< ": after update: " << track->get_cluster_IDs().back()
				<< endl;

				fout_chi2
				<< _event << "\t"
				<< iPHGenFitTrack << "\t"
				<< layer << "\t "
				<< iter->first
				<<endl;
#endif
			}
		}

		// Update other candidates
		if (incr_chi2s_new_tracks.size() > 1 and
				(layer >= _nlayers_maps and layer < _nlayers_maps + _nlayers_intt)) {
			for (auto iter = (++incr_chi2s_new_tracks.begin());
					iter != incr_chi2s_new_tracks.end(); ++iter) {

				if (!(iter->first < _max_splitting_chi2 and iter->first > 0))
					break;

#ifdef _DEBUG_
				std::cout <<__LINE__<<": "<< "Track Spliting with "<< "IncrChi2: "<< iter->first << std::endl;
#endif

//				_PHGenFitTracks.insert(
//						MapPHGenFitTrack::value_type(track_iter->first + iter->first,
//								std::shared_ptr < PHGenFit::Track> (iter->second)));

				_PHGenFitTracks.push_back(
						MapPHGenFitTrack::value_type(
								PHGenFitTrkProp::TrackQuality(
										tq.nhits + 1,
										tq.chi2  + iter->first,
										tq.ndf   + 2,
										tq.ntpc  + ((layer >= _nlayers_maps + _nlayers_intt) ? 1 : 0),
										tq.nintt + ((layer >= _nlayers_maps and layer < _nlayers_maps + _nlayers_intt) ? 1 : 0),
										tq.nmaps + ((layer < _nlayers_maps) ? 1 : 0)
										),
								std::shared_ptr < PHGenFit::Track> (iter->second)));
			}

#ifdef _DEBUG_
			std::cout <<__LINE__<<": "<< "_PHGenFitTracksSize: "<< _PHGenFitTracks.size() << std::endl;
			std::cout <<__LINE__<<": "<<track_iter->second->get_cluster_IDs().back() <<std::endl;
#endif
		}

#ifdef _DEBUG_
//		cout<<__LINE__<<": updateOneMeasurementKalman:"<<endl;
//		std::cout<<"iPHGenFitTrack: "<<iPHGenFitTrack
//				<<", layer: "<<layer
//				<<", #meas: "<<measurements.size()
//				<<", #tracks: "<<incr_chi2s_new_tracks.size()
//				<<", #totoal tracks: "<<_trackID_PHGenFitTrack.size()
//				<<std::endl;

		for (auto iter =
				incr_chi2s_new_tracks.begin();
				iter != incr_chi2s_new_tracks.end(); iter++) {
			std::cout << __LINE__ << ": IncrChi2: "<< iter->first << std::endl;
		}
#endif
		if(_analyzing_mode){
		  int ncand = 0;
		  for (auto iter =
			 incr_chi2s_new_tracks.begin();
		       iter != incr_chi2s_new_tracks.end(); iter++) {
		    if(iter->first<_max_incr_chi2s[layer] and iter->first > 0) ncand++;
		  }
		  /*
		    float this_pt = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Pt();
		    float this_phi = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Phi();
		    float this_eta = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom(state).Eta();
		  */
		  //"spt:seta:sphi:pt:eta:phi:layer:ncand:nmeas"
		  //_analyzing_ntuple->Fill(init_pt,init_eta,init_phi,this_pt,this_eta,this_phi,layer,ncand,measurements.size());
		}
		if(!layer_updated)
			++consecutive_missing_layer;
	} // layer loop

#ifdef _DEBUG_
	cout
	<< __LINE__
	<< ": clusterIDs size:  " << track->get_cluster_IDs().size()
	<< endl;
#endif

	//! Track succesfully propagated and return 0
	return 0;
}

PHGenFit::Measurement* PHGenFitTrkProp::SvtxClusterToPHGenFitMeasurement(
		const SvtxCluster* cluster) {

	if(!cluster) return nullptr;

	TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
	TVector3 n(cluster->get_x(), cluster->get_y(), 0);

	unsigned int begin_hit_id = *(cluster->begin_hits());

#ifdef _DEBUG_
	LogDebug(begin_hit_id)<<endl;
#endif
	SvtxHit* svtxhit = _svtxhitsmap->find(begin_hit_id)->second;

#ifdef _DEBUG_
	LogDebug(svtxhit->get_cellid())<<endl;
#endif
	PHG4Cell* cell_svtx = nullptr;
	PHG4Cell* cell_intt = nullptr;
	PHG4Cell* cell_maps = nullptr;

	if(_cells_svtx) cell_svtx = _cells_svtx->findCell(svtxhit->get_cellid());
	if(_cells_intt) cell_intt = _cells_intt->findCell(svtxhit->get_cellid());
	if(_cells_maps) cell_maps = _cells_maps->findCell(svtxhit->get_cellid());
	if(!(cell_svtx or cell_intt or cell_maps)){
		if(Verbosity()>=0)
			LogError("!(cell_svtx or cell_intt or cell_maps)")<<endl;
		return nullptr;
	}

	//17.4, 17.4, 17.4, 14.0, 14.0, 12.0, 11.5
	//float phi_tilt[7] = {0.304, 0.304, 0.304, 0.244, 0.244, 0.209, 0.201};

	unsigned int layer = cluster->get_layer();
	//std::cout << "cluster layer: " << layer << std::endl;
	if (cell_maps) {
		PHG4Cell* cell = cell_maps;

		int stave_index = cell->get_stave_index();
		int half_stave_index = cell->get_half_stave_index();
		int module_index = cell->get_module_index();
		int chip_index = cell->get_chip_index();

		double ladder_location[3] = { 0.0, 0.0, 0.0 };
		PHG4CylinderGeom_MVTX *geom =
				(PHG4CylinderGeom_MVTX*) _geom_container_maps->GetLayerGeom(
						layer);
		// returns the center of the sensor in world coordinates - used to get the ladder phi location
		geom->find_sensor_center(stave_index, half_stave_index,
				module_index, chip_index, ladder_location);
		//n.Print();
		n.SetXYZ(ladder_location[0], ladder_location[1], 0);
		n.RotateZ(geom->get_stave_phi_tilt());
		//n.Print();
	} else if (cell_intt) {
		PHG4Cell* cell = cell_intt;
		PHG4CylinderGeomINTT* geom =
		  dynamic_cast<PHG4CylinderGeomINTT*> (_geom_container_intt->GetLayerGeom(
							      layer));
		double hit_location[3] = { 0.0, 0.0, 0.0 };
		geom->find_segment_center(cell->get_ladder_z_index(),
				cell->get_ladder_phi_index(), hit_location);

		n.SetXYZ(hit_location[0], hit_location[1], 0);
		n.RotateZ(geom->get_strip_phi_tilt());
	}

	PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
			cluster->get_rphi_error(), cluster->get_z_error());

	meas->set_cluster_ID(cluster->get_id());

#ifdef _DEBUG_
	cout
	<<__LINE__
	<<": ID: " <<cluster->get_id()
	<<": layer: " << cluster->get_layer()
	<<": pos: {" << pos.X() <<", " << pos.Y() <<", " << pos.Z() << "}"
	<<": n: {" << n.X() <<", " << n.Y() <<", " <<n.Z() <<"}"
	//<<": rphi_error: " << cluster->get_rphi_error()
	<<": phi_error: " << cluster->get_phi_error()
	<<": theta error: " << cluster->get_z_error()/pos.Perp()
	<<endl;
#endif

	return meas;
}

int PHGenFitTrkProp::BuildLayerZPhiHitMap() {

	_layer_thetaID_phiID_cluserID.clear();

	//for(SimpleHit3D cluster : _clusters){
	for (SvtxClusterMap::Iter iter = _cluster_map->begin();
			iter != _cluster_map->end(); ++iter) {
		SvtxCluster* cluster = iter->second;

		unsigned int layer = cluster->get_layer();

		float x = cluster->get_x()-_vertex[0];
		float y = cluster->get_y()-_vertex[1];
		float z = cluster->get_z()-_vertex[2];

		float phi = atan2(y,x);
		float r = sqrt(x*x + y*y);
		float theta = atan2(r,z);

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

#ifdef _DEBUG_
//			cout
//			<<__LINE__<<": "
//			<<"{ "
//			<<layer <<", "
//			<<z <<", "
//			<<phi << "} =>"
//			<<idx << ": size: "
//			<<_layer_thetaID_phiID_cluserID.count(idx)
//			<<endl;
#endif

		_layer_thetaID_phiID_cluserID.insert(std::make_pair(idx, cluster->get_id()));
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<unsigned int> PHGenFitTrkProp::SearchHitsNearBy(const unsigned int layer,
		const float theta_center, const float phi_center, const float theta_window,
		const float phi_window) {

	std::vector<unsigned int> cluster_IDs;

	const unsigned int max_phi_bin = 16383; //2^14 - 1
	const unsigned int max_z_bin = 2047; // 2^11 - 1

	float lower_phi = phi_center-phi_window + _half_max_phi;
	float upper_phi = phi_center+phi_window + _half_max_phi;
	float lower_z = theta_center-theta_window + _half_max_theta;
	float upper_z = theta_center+theta_window + _half_max_theta;

	unsigned int lower_phi_bin = (unsigned int)((lower_phi)/_layer_thetaID_phiID_cluserID_phiSize);
	unsigned int upper_phi_bin = (unsigned int)((upper_phi)/_layer_thetaID_phiID_cluserID_phiSize);
	unsigned int lower_z_bin    = (unsigned int)(   (lower_z)/_layer_thetaID_phiID_cluserID_zSize);
	unsigned int upper_z_bin    = (unsigned int)(   (upper_z)/_layer_thetaID_phiID_cluserID_zSize);

	if(lower_phi < 0) lower_phi_bin= 0;
	if(upper_phi_bin > max_phi_bin) upper_phi_bin = max_phi_bin;

	if(lower_z < 0) lower_z_bin = 0;
	if(upper_z_bin > max_z_bin) upper_z_bin = max_z_bin;

	for (unsigned int iz = lower_z_bin; iz <= upper_z_bin; ++iz) {
		for (unsigned int irphi = lower_phi_bin; irphi <= upper_phi_bin;
				++irphi) {
			if(Verbosity() >= 2) _t_search_clusters_encoding->restart();
			unsigned int idx = encode_cluster_index(layer, iz, irphi);
			if(Verbosity() >= 2) _t_search_clusters_encoding->stop();

//#ifdef _DEBUG_
//			cout
//			<<__LINE__<<": "
//			<<"{ "
//			<<layer <<", "
//			<<iz <<", "
//			<<irphi << "} =>"
//			<<idx << ": size: "
//			<<_layer_thetaID_phiID_cluserID.count(idx)
//			<<endl;
//#endif
			if(Verbosity() >= 2) _t_search_clusters_map_iter->restart();
			for (auto iter = _layer_thetaID_phiID_cluserID.lower_bound(idx);
					iter != _layer_thetaID_phiID_cluserID.upper_bound(idx);
					++iter) {
				cluster_IDs.push_back(iter->second);
#ifdef _DEBUG_
				SvtxCluster* cluster = _cluster_map->get(iter->second);
				TVector3 v(cluster->get_x()-_vertex[0],cluster->get_y()-_vertex[1],cluster->get_z()-_vertex[2]);
				float phi_cluster = v.Phi();
				fout_kalman_pull
				<< _event << "\t"
				<< layer << "\t "
				<< phi_center - phi_cluster << "\t"
				<< theta_center - v.Theta() << "\t"
				<< phi_window/_search_wins_phi[layer] << "\t"
				<< theta_window/_search_wins_theta[layer] << "\t"
				<< (phi_center - phi_cluster)/phi_window*_search_wins_phi[layer] <<"\t "
				<< (theta_center - v.Theta())/theta_window*_search_wins_theta[layer] <<"\t"
				<< v.Perp()
				<<endl;
#endif
			}
			if(Verbosity() >= 2) _t_search_clusters_map_iter->stop();
		}
	}

#ifdef _DEBUG_
//	std::cout
//			<<__LINE__<<": "
//			<<"layer: "<<layer
//			<<": "<<phi_center <<" +- "<<phi_window
//			<<": "<<theta_center <<" +- "<<theta_window
//			<<endl;
	std::cout
			<<__LINE__<<": "
			<<"layer: "<<layer
			<<", rphi: {"<<lower_phi_bin
			<<", "<<upper_phi_bin
			<<"}, z: {"<<lower_z_bin
			<<", "<<upper_z_bin
			<<"}, found #clusters: "<<cluster_IDs.size()
			<<endl;
#endif

	return cluster_IDs;
}

unsigned int PHGenFitTrkProp::encode_cluster_index(const unsigned int layer,
		const float theta, const float phi) {

	unsigned int idx = UINT_MAX;

	if(layer >= 128) {
		LogError("layer >= 128\n");
		return idx;
	}

	if( (theta + _half_max_theta) < 0 ) {
		LogError("(theta + _half_max_theta) < 0 \n");
		return idx;
	}
	unsigned int itheta = (theta + _half_max_theta) / _layer_thetaID_phiID_cluserID_zSize;

	if( (phi + _half_max_phi) < 0 ) {
		LogError("(rphi + _half_max_rphi) < 0 \n");
		return idx;
	}
	unsigned int irphi = (phi + _half_max_phi) / _layer_thetaID_phiID_cluserID_phiSize;

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
		const unsigned int iz, const unsigned int irphi) {

	if(layer >= 128) {
		LogError("layer >= 128: ") << layer <<endl;;
		return UINT_MAX;
	}

	if(iz >= 2048) {
		LogError("iz >= 2048: ") << iz << endl;
		return UINT_MAX;
	}

	if(irphi >= 16384) {
		LogError("irphi >= 16384: ") << irphi <<endl;
		return UINT_MAX;
	}

	unsigned int index = 0;

	index |= (layer<<25);
	index |= (iz<<14);
	index |= (irphi);

	return index;
}

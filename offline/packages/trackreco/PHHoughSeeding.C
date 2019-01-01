/*!
 *  \file PHHoughSeeding.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#include "PHHoughSeeding.h"

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

PHHoughSeeding::PHHoughSeeding(
		const string& name,
		unsigned int nlayers_maps,
		unsigned int nlayers_intt,
		unsigned int nlayers_tpc,
		unsigned int nlayers_seeding,
		unsigned int min_nlayers_seeding
		)
    : PHTrackSeeding(name),
	  _t_seeding(nullptr),
	  _t_seed_init1(nullptr),
	  _t_seed_init2(nullptr),
	  _t_seed_init3(nullptr),
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
	  _seeding_layer(),
      _nlayers_seeding(nlayers_seeding),
      _min_nlayers_seeding(min_nlayers_seeding),
      _radii(),
      _material(),
      _user_material(),
      _magField(1.4),
      _reject_ghosts(true),
      _remove_hits(false),
      _min_pt(0.2),
      _min_z0(-14.0),
      _max_z0(+14.0),
      _max_r(1.0),
      _cut_on_dca(true),
      _dcaxy_cut(0.2),
      _dcaz_cut(0.2),
      _chi2_cut_fast_par0(10.0),
      _chi2_cut_fast_par1(50.0),
      _chi2_cut_fast_max(75.0),
      _chi2_cut_full(5.0),
      _ca_chi2_cut(5.0),
      _cos_angle_cut(0.985),
      _bin_scale(0.8),
      _z_bin_scale(0.8),
      _min_combo_hits(min_nlayers_seeding),
      _max_combo_hits(nlayers_seeding*4),
      _pt_rescale(0.9972 / 1.00117), // 1.0
      _fit_error_scale(_nlayers_seeding,1.0/sqrt(12.0)),
      _vote_error_scale(_nlayers_seeding,1.0),
      _layer_ilayer_map(),
      _clusters(),
      _tracks(),
      _track_errors(),
      _track_covars(),
      _vertex(),
      _tracker(NULL),
      _tracker_vertex(NULL),
      _tracker_etap_seed(NULL),
      _tracker_etam_seed(NULL),
      _vertexFinder(),
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
      _cut_min_pT(0.2),
      _nlayers_maps(nlayers_maps),
      _nlayers_intt(nlayers_intt),
      _nlayers_tpc(nlayers_tpc),
      _nlayers_all(_nlayers_maps+_nlayers_intt+_nlayers_tpc),
      _layer_ilayer_map_all(),
      _radii_all()
{
  _event = 0;
  
  _user_material.clear();
	for(unsigned int i=0;i<_nlayers_maps;++i)
		_user_material[i] = 0.003;
	for(unsigned int i=_nlayers_maps;i<_nlayers_maps+_nlayers_intt;++i)
		_user_material[i] = 0.008;

	//int seeding_layers[] = {7,15,25,35,45,55,66};
	int ninner_layer = _nlayers_maps+_nlayers_intt;
	int incr_layer = floor(_nlayers_tpc / 6.);
	int seeding_layers[] = {
			ninner_layer,
			ninner_layer+incr_layer*1,
			ninner_layer+incr_layer*2,
			ninner_layer+incr_layer*3,
			ninner_layer+incr_layer*4,
			ninner_layer+incr_layer*5,
			_nlayers_all-1};
	this->set_seeding_layer(seeding_layers, 7);

	_vertex_error.clear();
	_vertex_error.assign(3, 0.0100);

	if(_analyzing_mode){
	  cout << "Ana Mode, creating ntuples! " << endl;
	  _analyzing_file = new TFile("./PHHoughSeedingAna.root","RECREATE");
	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","pt:kappa:d:phi:dzdl:z0:nhit:ml:rec:dt");
	  cout << "Done" << endl;
	}
}

int PHHoughSeeding::Setup(PHCompositeNode* topNode) {

	// Start new interface ----
	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = PHTrackSeeding::Setup(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = CreateNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;
	// End new interface ----

	int min_layers    = 4;
	int nlayers_seeds = 7;
	int seeding_layers[] = {(int)(_nlayers_maps+_nlayers_intt),
				(int)(_nlayers_maps+_nlayers_intt+6),
				(int)(_nlayers_maps+_nlayers_intt+12),
				(int)(_nlayers_maps+_nlayers_intt+18),
				(int)(_nlayers_maps+_nlayers_intt+24),
				(int)(_nlayers_maps+_nlayers_intt+30),
				(int)(_nlayers_maps+_nlayers_intt+39)
				//7,13,19,25,31,37,46
	};
	
	set_seeding_layer(seeding_layers, nlayers_seeds);
	set_min_nlayers_seeding(min_layers);
	
	ret = InitializeGeometry(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK)
	  return ret;

	/*!
	 * Initilize parameters
	 */

	_t_seeding = new PHTimer("_t_seeding");
	_t_seeding->stop();

	_t_seed_init1 = new PHTimer("_t_seed_init1");
	_t_seed_init1->stop();

	_t_seed_init2 = new PHTimer("_t_seed_init2");
	_t_seed_init2->stop();

	_t_seed_init3 = new PHTimer("_t_seed_init3");
	_t_seed_init3->stop();

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

	if (Verbosity() > 0) {
		cout
				<< "====================== PHHoughSeeding::InitRun() ======================"
				<< endl;
		cout << " Magnetic field set to: " << _magField << " Tesla" << endl;
		cout << " Number of tracking layers: " << _nlayers_seeding << endl;
		for (unsigned int i = 0; i < _nlayers_seeding; ++i) {
			cout << "   Tracking layer #" << i << " " << "radius = "
					<< _radii[i] << " cm, " << "material = " << _material[i]
					<< endl;
			cout << "   Tracking layer #" << i << " " << "vote error scale = "
					<< _vote_error_scale[i] << ", " << "fit error scale = "
					<< _fit_error_scale[i] << endl;
		}
		cout << " Required hits: " << _min_nlayers_seeding << endl;
		cout << " Minimum pT: " << _min_pt << endl;
		cout << " Fast fit chisq cut min(par0+par1/pt,max): min( "
				<< _chi2_cut_fast_par0 << " + " << _chi2_cut_fast_par1
				<< " / pt, " << _chi2_cut_fast_max << " )" << endl;
		cout << " Maximum chisq (kalman fit): " << _chi2_cut_full << endl;
		cout << " Cell automaton chisq: " << _ca_chi2_cut << endl;
		cout << " Cos Angle Cut: " << _cos_angle_cut << endl;
		cout << " Ghost rejection: " << boolalpha << _reject_ghosts
				<< noboolalpha << endl;
		cout << " Hit removal: " << boolalpha << _remove_hits << noboolalpha
				<< endl;
		cout << " Maximum DCA: " << boolalpha << _cut_on_dca << noboolalpha
				<< endl;
		if (_cut_on_dca) {
			cout << "   Maximum DCA cut: " << _dcaxy_cut << endl;
		}
		cout << "   Maximum DCAZ cut: " << _dcaz_cut << endl;
		cout << " Phi bin scale: " << _bin_scale << endl;
		cout << " Z bin scale: " << _z_bin_scale << endl;
		cout << " Momentum rescale factor: " << _pt_rescale << endl;
		cout
				<< "==========================================================================="
				<< endl;
	}
	
	return ret;
}

int PHHoughSeeding::Process() {

  if (Verbosity() > 0){
	  cout << "PHHoughSeeding::process_event -- entered" << endl;
	  cout << "nMapsLayers = " << _nlayers_maps << endl;
	  cout << "nInttLayers = " << _nlayers_intt << endl;
	  cout << "nTPCLayers = " << _nlayers_tpc << endl;
  }
	// start fresh
	int code;

	_clusters.clear();
	_all_tracks.clear();
	_all_track_errors.clear();
	_all_track_covars.clear();

	_vertex.clear();
	_vertex.assign(3, 0.0);
	
	/*iteration from Christof*/ {
	  _tracks.clear();
	  _track_errors.clear();
	  _track_covars.clear();
	  
	  /*iteration 1 from Christof*/ {
	    int min_layers    = 4;
	    int nlayers_seeds = 7;
	    int seeding_layers[] = {(int)(_nlayers_maps+_nlayers_intt),
				    (int)(_nlayers_maps+_nlayers_intt+8),
				    (int)(_nlayers_maps+_nlayers_intt+16),
				    (int)(_nlayers_maps+_nlayers_intt+24),
				    (int)(_nlayers_maps+_nlayers_intt+32),
				    (int)(_nlayers_maps+_nlayers_intt+40),
				    (int)(_nlayers_maps+_nlayers_intt+45)  // avoid the outer TPC layer, it is inefficient
				    //7,13,19,25,31,37,46
	    };

	    set_seeding_layer(seeding_layers, nlayers_seeds);
	    set_min_nlayers_seeding(min_layers);
	    _min_combo_hits = min_layers;
	    _max_combo_hits = nlayers_seeds;

	    if(Verbosity() >= 1) _t_seed_init1->restart();
	  }

	  _min_nlayers_seeding--;
	  if(Verbosity() >= 1) _t_seeding->restart();
	  
	  //-----------------------------------
	  // Translate into Helix_Hough objects
	  //-----------------------------------
	  
	  code = translate_input();//Check if cluster is already used in here
	  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

	  //-----------------------------------
	  // Obtain initial vertex from SvtxVtxMap[0]
	  //-----------------------------------
	  {
	    code = vertexing();
	    if (code != Fun4AllReturnCodes::EVENT_OK) return code;
	    // here expect vertex to be better than +/- 500 um
	  }

	  //-----------------------------------
	  // Seeding - Alan's Hough Tracking with selected layers
	  //-----------------------------------
	  code = full_track_seeding();
	  if (code != Fun4AllReturnCodes::EVENT_OK)
	    return code;
	  
	  if(Verbosity() >= 1) _t_seeding->stop();
	  
	  _t_seed_init1->stop();

	  if(Verbosity() > 1) print_timers();

	} /*end of iteration*/

	//	CleanupTracksByHitPattern();

  //-----------------------------------
  // Alan's exportation
  //-----------------------------------

	code = export_output();
	if (code != Fun4AllReturnCodes::EVENT_OK)
	  return code;

	++_event;
	
	return Fun4AllReturnCodes::EVENT_OK;
}


void PHHoughSeeding::print_timers() {
  
  std::cout << "=============== Timers: ===============" << std::endl;
  std::cout << "CPUSCALE Seeding time:                "<<_t_seeding->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "CPUSCALE Init Seed1 time:                "<<_t_seed_init1->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "CPUSCALE Init Seed2 time:                "<<_t_seed_init2->get_accumulated_time()/1000. << " sec" <<std::endl;
  std::cout << "CPUSCALE Init Seed3 time:                "<<_t_seed_init3->get_accumulated_time()/1000. << " sec" <<std::endl;
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

int PHHoughSeeding::End() {

	delete _t_seeding;
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

	delete _tracker_etap_seed;
	_tracker_etap_seed = NULL;
	delete _tracker_etam_seed;
	_tracker_etam_seed = NULL;
	delete _tracker_vertex;
	_tracker_vertex = NULL;
	delete _tracker;
	_tracker = NULL;


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


void PHHoughSeeding::set_material(int layer, float value) {
	_user_material[layer] = value;
}

int PHHoughSeeding::CreateNodes(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::InitializeGeometry(PHCompositeNode *topNode) {
  
  //---------------------------------------------------------
  // Grab Run-Dependent Detector Geometry and Configure Hough
  //---------------------------------------------------------
  
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<
  PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  
  _nlayers_seeding = _seeding_layer.size();
	
  //=================================================//
  //  Initializing HelixHough objects                //
  //=================================================//
  
  // Since the G4 layers don't necessarily correspond to the
  // silicon layers, and don't necessarily start from zero (argh),
  // we create our own layers numbers that are consecutive
  // starting from zero.
  
  // Now that we have two kinds of layers, I won't know in principle
  // which type is in what order, so I figure that out now...

  _radii.assign(_nlayers_seeding, 0.0);
  map<float, int> radius_layer_map;
  
  _radii_all.assign(_nlayers_all, 0.0);
  _layer_ilayer_map.clear();
  _layer_ilayer_map_all.clear();
  if (cellgeos) {
    PHG4CylinderCellGeomContainer::ConstRange layerrange =
      cellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
	   layerrange.first; layeriter != layerrange.second; ++layeriter) {
      radius_layer_map.insert(
			      make_pair(layeriter->second->get_radius(),
					layeriter->second->get_layer()));
    }
  }
  
  if (laddergeos) {
    PHG4CylinderGeomContainer::ConstRange layerrange =
      laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
	   layerrange.first; layeriter != layerrange.second; ++layeriter) {
      radius_layer_map.insert(
			      make_pair(layeriter->second->get_radius(),
					layeriter->second->get_layer()));
    }
  }
  
  if (mapsladdergeos) {
    PHG4CylinderGeomContainer::ConstRange layerrange =
      mapsladdergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
	   layerrange.first; layeriter != layerrange.second; ++layeriter) {
      radius_layer_map.insert(
			      make_pair(layeriter->second->get_radius(),
					layeriter->second->get_layer()));
    }
  }
  
  //	if (Verbosity() >= 2) {
  //		for (map<float, int>::const_iterator iter = radius_layer_map.begin();
  //				iter != radius_layer_map.end(); iter++) {
  //			cout << "radius_layer_map: first: " << iter->first << "; second: "
  //					<< iter->second << endl;
  //		}
  //	}
  
  // now that the layer ids are sorted by radius, I can create a storage
  // index, ilayer, that is 0..N-1 and sorted by radius
  
  int ilayer = 0;
  for (map<float, int>::iterator iter = radius_layer_map.begin();
       iter != radius_layer_map.end(); ++iter) {
    _layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));
    
    if (std::find(_seeding_layer.begin(), _seeding_layer.end(),
		  iter->second) != _seeding_layer.end()) {
      _layer_ilayer_map.insert(make_pair(iter->second, ilayer));
      ++ilayer;
    }
    //if(ilayer >= (int) _radii.size()) break; //yuhw
  }
  
  //	if (Verbosity() >= 10) {
  //		for (map<int, unsigned int>::const_iterator iter = _layer_ilayer_map_all.begin();
  //				iter != _layer_ilayer_map_all.end(); iter++) {
  //			cout << "_layer_ilayer_map_all: first: " << iter->first << "; second: "
  //					<< iter->second << endl;
  //		}
  //	}
  
  // now we extract the information from the cellgeos first
  if (cellgeos) {
    PHG4CylinderCellGeomContainer::ConstRange begin_end =
      cellgeos->get_begin_end();
    PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++) {
      PHG4CylinderCellGeom *geo = miter->second;
      
      //if(cellgeo->get_layer() > (int) _radii.size() ) continue;
      
      //			if (Verbosity() >= 2)
      //				cellgeo->identify();
      
      //TODO
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
	geo->get_radius() + 0.5 * geo->get_thickness();
      
      
      if (_layer_ilayer_map.find(geo->get_layer())
	  != _layer_ilayer_map.end()) {
	_radii[_layer_ilayer_map[geo->get_layer()]] =
	  geo->get_radius();
      }
    }
  }
  
  if (laddergeos) {
    PHG4CylinderGeomContainer::ConstRange begin_end =
      laddergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++) {
      PHG4CylinderGeom *geo = miter->second;
      
      //if(geo->get_layer() > (int) _radii.size() ) continue;
      
      //			if (Verbosity() >= 2)
      //				geo->identify();
      
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
	geo->get_radius() + 0.5*geo->get_thickness();
      
      if (_layer_ilayer_map.find(geo->get_layer())
	  != _layer_ilayer_map.end()) {
	_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
      }
    }
  }
  
  if (mapsladdergeos) {
    PHG4CylinderGeomContainer::ConstRange begin_end =
      mapsladdergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; miter++) {
      PHG4CylinderGeom *geo = miter->second;
      
      //if(geo->get_layer() > (int) _radii.size() ) continue;
      
      //			if (Verbosity() >= 2)
      //				geo->identify();
      
      //TODO
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
	geo->get_radius();
      
      if (_layer_ilayer_map.find(geo->get_layer())
	  != _layer_ilayer_map.end()) {
	_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
      }
    }
  }
  // set material on each layer
  
  _material.assign(_radii.size(), 0.03);
  
  map<int, float>::iterator mat_it;
  for (map<int, float>::iterator iter = _user_material.begin();
       iter != _user_material.end(); ++iter) {
    if (_layer_ilayer_map.find(iter->first) != _layer_ilayer_map.end()) {
      _material[_layer_ilayer_map[iter->first]] = iter->second;
    }
  }
  if(_tracker) delete _tracker;
  if(_tracker_vertex) delete _tracker_vertex;
  if(_tracker_etap_seed) delete _tracker_etap_seed;
  
  // initialize the pattern recogition tools
  setup_tracker_object();
  setup_initial_tracker_object();
  setup_seed_tracker_objects();
  
  
  /*!
   * Now have to load geometry nodes to get norm vector
   */
  
  _cells_svtx = findNode::getClass<PHG4CellContainer>(topNode,
						      "G4CELL_SVTX");
  
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

int PHHoughSeeding::setup_seed_tracker_objects() {

	float kappa_max = ptToKappa(_min_pt);

	// for the initial tracker we may not have the best guess on the vertex yet
	// so I've doubled the search range on dca and dcaz

	std::vector<unsigned int> onezoom(5, 0);
	std::vector<vector<unsigned int> > zoomprofile;
	zoomprofile.assign(5, onezoom);
	zoomprofile[0][0] = 16;
	zoomprofile[0][1] = 1;
	zoomprofile[0][2] = 4;
	zoomprofile[0][3] = 8;
	zoomprofile[0][4] = 1;

	zoomprofile[1][0] = 16;
	zoomprofile[1][1] = 1;
	zoomprofile[1][2] = 4;
	zoomprofile[1][3] = 4;
	zoomprofile[1][4] = 2;

	zoomprofile[2][0] = 4;
	zoomprofile[2][1] = 3;
	zoomprofile[2][2] = 2;
	zoomprofile[2][3] = 1;
	zoomprofile[2][4] = 3;

	for (unsigned int i = 2; i <= 3; ++i) {
		zoomprofile[i][0] = 3;
		zoomprofile[i][1] = 3;
		zoomprofile[i][2] = 3;
		zoomprofile[i][3] = 3;
		zoomprofile[i][4] = 3;
	}

	HelixRange pos_range(0.0, 2. * M_PI,  // center of rotation azimuthal angles
			-_max_r, _max_r,   // 2d dca range
			0.0, kappa_max,    // curvature range
			0.0, 0.9,          // dzdl range
			_min_z0, _max_z0); // dca_z range

	_tracker_etap_seed = new sPHENIXSeedFinder(zoomprofile, 1, pos_range,
			_material, _radii, _magField);
	_tracker_etap_seed->setNLayers(_nlayers_seeding);
	_tracker_etap_seed->requireLayers(_min_nlayers_seeding);
	_tracker_etap_seed->setClusterStartBin(1);
	_tracker_etap_seed->setRejectGhosts(_reject_ghosts);
	_tracker_etap_seed->setFastChi2Cut(_chi2_cut_fast_par0, _chi2_cut_fast_par1,
			_chi2_cut_fast_max);
	_tracker_etap_seed->setChi2Cut(_chi2_cut_full);
	_tracker_etap_seed->setChi2RemovalCut(_chi2_cut_full * 0.5);
	_tracker_etap_seed->setCellularAutomatonChi2Cut(_ca_chi2_cut);
	_tracker_etap_seed->setPrintTimings(false);
	_tracker_etap_seed->setCutOnDca(false);
	_tracker_etap_seed->setSmoothBack(true);
	_tracker_etap_seed->setBinScale(_bin_scale);
	_tracker_etap_seed->setZBinScale(_z_bin_scale);
	_tracker_etap_seed->setRemoveHits(_remove_hits);
	_tracker_etap_seed->setSeparateByHelicity(true);
	_tracker_etap_seed->setMaxHitsPairs(0);
	_tracker_etap_seed->setCosAngleCut(_cos_angle_cut);

	for (unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
		float scale1 = _fit_error_scale[ilayer];
		float scale2 = _vote_error_scale[ilayer];
		float scale = scale1 / scale2;
		_tracker_etap_seed->setHitErrorScale(ilayer, scale);
	}

	// for the initial tracker we may not have the best guess on the vertex yet
	// so I've doubled the search range on dca and dcaz

	HelixRange neg_range(0.0, 2. * M_PI,  // center of rotation azimuthal angles
			-_max_r, _max_r,   // 2d dca range
			0.0, kappa_max,    // curvature range
			-0.9, 0.0,         // dzdl range
			_min_z0, _max_z0); // dca_z range

	_tracker_etam_seed = new sPHENIXSeedFinder(zoomprofile, 1, neg_range,
			_material, _radii, _magField);
	_tracker_etam_seed->setNLayers(_nlayers_seeding);
	_tracker_etam_seed->requireLayers(_min_nlayers_seeding);
	_tracker_etam_seed->setClusterStartBin(1);
	_tracker_etam_seed->setRejectGhosts(_reject_ghosts);
	_tracker_etam_seed->setFastChi2Cut(_chi2_cut_fast_par0, _chi2_cut_fast_par1,
			_chi2_cut_fast_max);
	_tracker_etam_seed->setChi2Cut(_chi2_cut_full);
	_tracker_etam_seed->setChi2RemovalCut(_chi2_cut_full * 0.5);
	_tracker_etam_seed->setCellularAutomatonChi2Cut(_ca_chi2_cut);
	_tracker_etam_seed->setPrintTimings(false);
	_tracker_etam_seed->setCutOnDca(false);
	_tracker_etam_seed->setSmoothBack(true);
	_tracker_etam_seed->setBinScale(_bin_scale);
	_tracker_etam_seed->setZBinScale(_z_bin_scale);
	_tracker_etam_seed->setRemoveHits(_remove_hits);
	_tracker_etam_seed->setSeparateByHelicity(true);
	_tracker_etam_seed->setMaxHitsPairs(0);
	_tracker_etam_seed->setCosAngleCut(_cos_angle_cut);

	for (unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
		float scale1 = _fit_error_scale[ilayer];
		float scale2 = _vote_error_scale[ilayer];
		float scale = scale1 / scale2;
		_tracker_etam_seed->setHitErrorScale(ilayer, scale);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::setup_initial_tracker_object() {

	// copy of the final tracker modified to:
	// expand the DCA search regions (2.0 cm z search > 3 sigma of BBC z vertex
	// remove the DCA cut on the track output

	// tell the initial tracker object the phase space extent of the search region
	// and the recursive zoom factors to utilize

	float kappa_max = ptToKappa(_min_pt);

	// for the initial tracker we may not have the best guess on the vertex yet
	// so I've doubled the search range on dca and dcaz

	HelixRange top_range(0.0, 2. * M_PI,  // center of rotation azimuthal angles
			-1.0, +1.0,      // 2d dca range
			0.0, kappa_max,  // curvature range
			-0.9, 0.9,       // dzdl range
			-2.0, +2.0);     // dca_z range

	vector<unsigned int> onezoom(5, 0);
	vector<vector<unsigned int> > zoomprofile;
	zoomprofile.assign(5, onezoom);
	zoomprofile[0][0] = 16;
	zoomprofile[0][1] = 1;
	zoomprofile[0][2] = 4;
	zoomprofile[0][3] = 8;
	zoomprofile[0][4] = 1;

	zoomprofile[1][0] = 16;
	zoomprofile[1][1] = 1;
	zoomprofile[1][2] = 4;
	zoomprofile[1][3] = 4;
	zoomprofile[1][4] = 2;

	zoomprofile[2][0] = 4;
	zoomprofile[2][1] = 3;
	zoomprofile[2][2] = 2;
	zoomprofile[2][3] = 1;
	zoomprofile[2][4] = 3;

	for (unsigned int i = 2; i <= 3; ++i) {
		zoomprofile[i][0] = 3;
		zoomprofile[i][1] = 3;
		zoomprofile[i][2] = 3;
		zoomprofile[i][3] = 3;
		zoomprofile[i][4] = 3;
	}

	_tracker_vertex = new sPHENIXSeedFinder(zoomprofile, 1, top_range,
			_material, _radii, _magField);
	_tracker_vertex->setNLayers(_nlayers_seeding);
	_tracker_vertex->requireLayers(_min_nlayers_seeding);
	_tracker_vertex->setClusterStartBin(1);
	_tracker_vertex->setRejectGhosts(_reject_ghosts);
	_tracker_vertex->setFastChi2Cut(_chi2_cut_fast_par0, _chi2_cut_fast_par1,
			_chi2_cut_fast_max);
	_tracker_vertex->setChi2Cut(_chi2_cut_full);
	_tracker_vertex->setChi2RemovalCut(_chi2_cut_full * 0.5);
	_tracker_vertex->setCellularAutomatonChi2Cut(_ca_chi2_cut);
	_tracker_vertex->setPrintTimings(false);
	_tracker_vertex->setCutOnDca(false);
	_tracker_vertex->setSmoothBack(true);
	_tracker_vertex->setBinScale(_bin_scale);
	_tracker_vertex->setZBinScale(_z_bin_scale);
	_tracker_vertex->setRemoveHits(_remove_hits);
	_tracker_vertex->setSeparateByHelicity(true);
	_tracker_vertex->setMaxHitsPairs(0);
	_tracker_vertex->setCosAngleCut(_cos_angle_cut);

	for (unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
		float scale1 = _fit_error_scale[ilayer];
		float scale2 = _vote_error_scale[ilayer];
		float scale = scale1 / scale2;
		_tracker_vertex->setHitErrorScale(ilayer, scale);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::setup_tracker_object() {

	// input vertex must be within 500 um of final

	// tell the tracker object the phase space extent of the search region
	// and the recursive zoom factors to utilize

	float kappa_max = ptToKappa(_min_pt);

	HelixRange top_range(0.0, 2. * M_PI,  // center of rotation azimuthal angles
			-_dcaxy_cut, _dcaxy_cut,  // 2d dca range
			0.0, kappa_max,           // curvature range
			-0.9, 0.9,                // dzdl range
			-_dcaz_cut, _dcaz_cut);   // dca_z range

	vector<unsigned int> onezoom(5, 0);
	vector<vector<unsigned int> > zoomprofile;
	zoomprofile.assign(5, onezoom);
	zoomprofile[0][0] = 16;
	zoomprofile[0][1] = 1;
	zoomprofile[0][2] = 4;
	zoomprofile[0][3] = 8;
	zoomprofile[0][4] = 1;

	zoomprofile[1][0] = 16;
	zoomprofile[1][1] = 1;
	zoomprofile[1][2] = 4;
	zoomprofile[1][3] = 4;
	zoomprofile[1][4] = 2;

	zoomprofile[2][0] = 4;
	zoomprofile[2][1] = 3;
	zoomprofile[2][2] = 2;
	zoomprofile[2][3] = 1;
	zoomprofile[2][4] = 3;

	for (unsigned int i = 2; i <= 3; ++i) {
		zoomprofile[i][0] = 3;
		zoomprofile[i][1] = 3;
		zoomprofile[i][2] = 3;
		zoomprofile[i][3] = 3;
		zoomprofile[i][4] = 3;
	}

	_tracker = new sPHENIXSeedFinder(zoomprofile, 1, top_range, _material,
			_radii, _magField);
	_tracker->setNLayers(_nlayers_seeding);
	_tracker->requireLayers(_min_nlayers_seeding);
	_tracker->setClusterStartBin(1);
	_tracker->setRejectGhosts(_reject_ghosts);
	_tracker->setFastChi2Cut(_chi2_cut_fast_par0, _chi2_cut_fast_par1,
			_chi2_cut_fast_max);
	_tracker->setChi2Cut(_chi2_cut_full);
	_tracker->setChi2RemovalCut(_chi2_cut_full * 0.5);
	_tracker->setCellularAutomatonChi2Cut(_ca_chi2_cut);
	_tracker->setPrintTimings(false);
	if(Verbosity() >= 2)
		_tracker->setPrintTimings(true);
	//_tracker->setVerbosity(Verbosity());
	_tracker->setCutOnDca(_cut_on_dca);
	_tracker->setDcaCut(_dcaxy_cut);
	_tracker->setSmoothBack(true);
	_tracker->setBinScale(_bin_scale);
	_tracker->setZBinScale(_z_bin_scale);
	_tracker->setRemoveHits(_remove_hits);
	_tracker->setSeparateByHelicity(true);
	_tracker->setMaxHitsPairs(0);
	_tracker->setCosAngleCut(_cos_angle_cut);

	for (unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
		float scale1 = _fit_error_scale[ilayer];
		float scale2 = _vote_error_scale[ilayer];
		float scale = scale1 / scale2;
		_tracker->setHitErrorScale(ilayer, scale);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::GetNodes(PHCompositeNode* topNode) {

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

//	if(_hit_used_map_size!=0) delete[] _hit_used_map;
//	_hit_used_map_size = static_cast<int>(_cluster_map->size());
//	_hit_used_map = new int[_hit_used_map_size];
//	for (Int_t i=0;i<_hit_used_map_size;i++){
//	  _hit_used_map[i] = 0;
//	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::translate_input() {

  _clusters.clear();
  int count = 0;
  int count7 = 0;
  int count46 = 0;
  int nhits[60];
  int nhits_all[60];
  for(int i = 0; i< 60 ;i++){
     nhits[i] = 0;
     nhits_all[i] = 0;
  }
  for (SvtxClusterMap::Iter iter = _cluster_map->begin();
       iter != _cluster_map->end(); ++iter) {
    //if(_hit_used_map[iter->first]!=0){continue;}
    if(_assoc_container->GetTracksFromCluster(iter->first).size()>0){continue;}
    count++;
    SvtxCluster* cluster = iter->second;
    nhits_all[cluster->get_layer()]++;
    if(cluster->get_layer()==(unsigned int)(_nlayers_maps+_nlayers_intt))count7++;
    if(cluster->get_layer()==(unsigned int)(_nlayers_maps+_nlayers_intt+40))count46++;
    //	  cout << "first: " << iter->first << endl; 
    /*
      float vz = 0.0;
      float x  = cluster->get_x();
      float y  = cluster->get_y();
      float z  = cluster->get_z();
      float dz = z - vz;
      float r  = sqrt(x*x+y*y);
      float zsize = cluster->get_z_size();
      bool goodhit = false;
      
      if(TMath::Abs(dz)<40&&zsize<3)
      goodhit = true;
      
      if(zsize > (TMath::Abs(dz)/r * 2.448 + 0.5)&&
      zsize < (TMath::Abs(dz)/r * 2.448 + 3.5) )
      goodhit = true;
      if(goodhit==false) continue;
      //ntp_cluster.Draw("zsize:z-gvz","layer==7&&zsize>(abs(z-gvz)*0.08+0.5)&&zsize<(abs(z-gvz)*0.08)+3.5")
      */
    //unsigned int ilayer = _layer_ilayer_map[cluster->get_layer()];
    
    //		unsigned int ilayer = _layer_ilayer_map_all[cluster->get_layer()];
    //		if(ilayer >= _nlayers_seeding) continue;
    
    unsigned int ilayer = UINT_MAX;
    std::map<int, unsigned int>::const_iterator it = _layer_ilayer_map.find(cluster->get_layer());
    if(it != _layer_ilayer_map.end())
      ilayer = it->second;
    if(ilayer >= _nlayers_seeding) continue;
    
    SimpleHit3D hit3d;
    
    hit3d.set_id(cluster->get_id());
    hit3d.set_layer(ilayer);
    
    hit3d.set_x(cluster->get_x());
    hit3d.set_y(cluster->get_y());
    hit3d.set_z(cluster->get_z());
    
    // hit3d.set_ex(2.0*sqrt(cluster->get_size(0,0)));
    // hit3d.set_ey(2.0*sqrt(cluster->get_size(1,1)));
    // hit3d.set_ez(2.0*sqrt(cluster->get_size(2,2)));
    
    // copy covariance over
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
	hit3d.set_error(i, j, cluster->get_error(i, j));
	
	//FIXME
	//hit3d.set_size(i, j, cluster->get_size(i, j)); // original
	hit3d.set_size(i, j, cluster->get_error(i, j)*sqrt(12.)); // yuhw 2017-05-08
      }
    }
    /*    float x  = cluster->get_x();
    float y  = cluster->get_y();
    float z  = cluster->get_z();
    float r  = sqrt(x*x+y*y);
    */
    nhits[ilayer]++;
    _clusters.push_back(hit3d);
  }
  
  if (Verbosity() > 20) {
    cout
      << "-------------------------------------------------------------------"
      << endl;
    cout
      << "PHHoughSeeding::process_event has the following input clusters:"
      << endl;
    
    for (unsigned int i = 0; i < _clusters.size(); ++i) {
      cout << "n init clusters = " << _clusters.size() << endl;
      _clusters[i].print();
    }
    
    cout
      << "-------------------------------------------------------------------"
      << endl;
  }
  
  if(Verbosity() >= 1){
    cout << "CPUSCALE hits: " << count << endl;
    }
  if(Verbosity() >= 10){
    for(int i  = 0;i<60;i++){
      cout << "layer: " << i << " << hits: " << nhits[i] << " | " << nhits_all[i] << endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::fast_vertex_from_bbc() {

	// fail over to bbc vertex if no tracks were found...
	if (_bbc_vertexes) {

		BbcVertex* vertex = _bbc_vertexes->begin()->second;

		if (vertex) {

			_vertex[0] = 0.0;
			_vertex[1] = 0.0;
			_vertex[2] = vertex->get_z();

			if (Verbosity())
				cout << " initial bbc vertex guess: " << _vertex[0] << " "
						<< _vertex[1] << " " << _vertex[2] << endl;
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::fast_vertex_guessing() {

	// fast vertex guessing uses two tracker objects
	// one looks for postive going eta tracks, the other for negative going tracks
	// it asks for just a handful of each, but searches
	// over the broadest possible phase space for the vertex origin

	// limit each window to no more than N tracks
	unsigned int maxtracks = 10;

	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	// loop over initial tracking windows
	std::vector<SimpleTrack3D> newtracks;

	_tracker_etap_seed->clear();
	_tracker_etap_seed->findHelices(_clusters, _min_combo_hits, _max_combo_hits,
			newtracks, maxtracks);

	for (unsigned int t = 0; t < newtracks.size(); ++t) {
		_tracks.push_back(newtracks[t]);
		_track_covars.push_back((_tracker_etap_seed->getKalmanStates())[t].C);
	}

	_tracker_etap_seed->clear();
	newtracks.clear();

	_tracker_etam_seed->clear();
	_tracker_etam_seed->findHelices(_clusters, _min_combo_hits, _max_combo_hits,
			newtracks, maxtracks);

	for (unsigned int t = 0; t < newtracks.size(); ++t) {
		_tracks.push_back(newtracks[t]);
		_track_covars.push_back((_tracker_etam_seed->getKalmanStates())[t].C);
	}

	_tracker_etam_seed->clear();
	newtracks.clear();

	_vertex.clear();
	_vertex.assign(3, 0.0);

	if (Verbosity())
		cout << " seed track finding count: " << _tracks.size() << endl;

	if (!_tracks.empty()) {

		// --- compute seed vertex from initial track finding --------------------

		double zsum = 0.0;

		for (unsigned int i = 0; i < _tracks.size(); ++i) {
			zsum += _tracks[i].z0;
		}

		_vertex[2] = zsum / _tracks.size();

		if (Verbosity() > 0) {
			cout << " seed track vertex pre-fit: " << _vertex[0] << " "
					<< _vertex[1] << " " << _vertex[2] << endl;
		}

		// start with the average position and converge from there
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 3.00, true);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.10, true);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.02, true);
	}

	// we don't need the tracks anymore
	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	if (Verbosity() > 0) {
		cout << " seed track vertex post-fit: " << _vertex[0] << " "
				<< _vertex[1] << " " << _vertex[2] << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::initial_vertex_finding() {

	// shift to the guess vertex position
	// run the tracking pattern recognition, stop after some number of tracks
	// have been found, fit those tracks to a vertex
	// nuke out tracks, leave vertex info, shift back

	float shift_dx = -_vertex[0];
	float shift_dy = -_vertex[1];
	float shift_dz = -_vertex[2];

	// shift to vertex guess position
	shift_coordinate_system(shift_dx, shift_dy, shift_dz);

	// limit each window to no more than N tracks
	unsigned int maxtracks = 40;

	// reset track storage and tracker
	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	_tracker_vertex->clear();

	// initial track finding
	_tracker_vertex->findHelices(_clusters, _min_combo_hits, _max_combo_hits,
			_tracks, maxtracks);

	for (unsigned int t = 0; t < _tracks.size(); ++t) {
		_track_covars.push_back((_tracker_vertex->getKalmanStates())[t].C);
	}

	// don't need the tracker object anymore
	_tracker_vertex->clear();

	if (Verbosity())
		cout << " initial track finding count: " << _tracks.size() << endl;

	if (!_tracks.empty()) {

		// --- compute seed vertex from initial track finding --------------------

		double xsum = 0.0;
		double ysum = 0.0;
		double zsum = 0.0;

		for (unsigned int i = 0; i < _tracks.size(); ++i) {
			xsum += _tracks[i].d * cos(_tracks[i].phi);
			ysum += _tracks[i].d * sin(_tracks[i].phi);
			zsum += _tracks[i].z0;
		}

		_vertex[0] = xsum / _tracks.size();
		_vertex[1] = ysum / _tracks.size();
		_vertex[2] = zsum / _tracks.size();

		if (Verbosity() > 0) {
			cout << " initial track vertex pre-fit: " << _vertex[0] - shift_dx
					<< " " << _vertex[1] - shift_dy << " "
					<< _vertex[2] - shift_dz << endl;
		}

		// start with the average position and converge from there
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 3.00, true);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.10, true);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.02, false);

	}

	// don't need the tracks anymore
	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	// shift back to the global coordinates
	shift_coordinate_system(-shift_dx, -shift_dy, -shift_dz);

	if (Verbosity() > 0) {
		cout << " initial track vertex post-fit: " << _vertex[0] << " "
				<< _vertex[1] << " " << _vertex[2] << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}


int PHHoughSeeding::vertexing() {

	SvtxVertex * svtx_vtx = _vertex_map->get(0);

	_vertex.clear();
	_vertex.assign(3, 0.0);

	_vertex[0] = svtx_vtx->get_x();
	_vertex[1] = svtx_vtx->get_y();
	_vertex[2] = svtx_vtx->get_z();

	_vertex_error[0] = sqrt(svtx_vtx->get_error(0,0));
	_vertex_error[1] = sqrt(svtx_vtx->get_error(1,2));
	_vertex_error[2] = sqrt(svtx_vtx->get_error(2,2));


	if (Verbosity() > 1) {
		cout << __LINE__ << " PHHoughSeeding::vertexing: {" << _vertex[0]
				<< ", " << _vertex[1] << ", " << _vertex[2] << "} +- {"
				<< _vertex_error[0] << ", " << _vertex_error[1] << ", "
				<< _vertex_error[2] << "}" << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::full_track_seeding() {

	float shift_dx = -_vertex[0];
	float shift_dy = -_vertex[1];
	float shift_dz = -_vertex[2];
	// shift to initial vertex position
	shift_coordinate_system(shift_dx, shift_dy, shift_dz);

	// reset track storage and tracker
	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	_tracker->clear();
	// final track finding
	_tracker->findHelices(_clusters, _min_combo_hits, _max_combo_hits, _tracks);
	if(Verbosity() >= 1)
	  cout << "SEEDSTUDY nbefore clean (" << _min_nlayers_seeding << "): " << _tracks.size() << endl;
	// Cleanup Seeds
#ifdef _USE_ALAN_TRACK_REFITTING_
#else
	if(Verbosity() >= 1) _t_seeds_cleanup->restart();
	CleanupSeeds();
	if(Verbosity() >= 1) _t_seeds_cleanup->stop();
#endif
	if(Verbosity() >= 1)
	  cout << "SEEDSTUDY nafter clean: " << _tracks.size() << endl;
	for (unsigned int tt = 0; tt < _tracks.size(); ++tt) {
		_track_covars.push_back((_tracker->getKalmanStates())[tt].C);
		_track_errors.push_back(_tracker->getKalmanStates()[tt].chi2);
	}

	// we will need the tracker object below to refit the tracks... so we won't
	// reset it here

	if (Verbosity() > 0)
		cout << " final track count: " << _tracks.size() << endl;
#ifdef _USE_ALAN_FULL_VERTEXING_
	if (!_tracks.empty()) {

		if (Verbosity() > 0) {
			cout << " final vertex pre-fit: " << _vertex[0] - shift_dx << " "
					<< _vertex[1] - shift_dy << " " << _vertex[2] - shift_dz
					<< endl;
		}

		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.300, false);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.100, false);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.020, false);
		_vertexFinder.findVertex(_tracks, _track_covars, _vertex, 0.005, false);

		if (Verbosity() > 0) {
			cout << " final vertex post-fit: " << _vertex[0] - shift_dx << " "
					<< _vertex[1] - shift_dy << " " << _vertex[2] - shift_dz
					<< endl;
		}
	}
#endif
	// shift back to global coordinates
	shift_coordinate_system(-shift_dx, -shift_dy, -shift_dz);

#ifdef _USE_ALAN_TRACK_REFITTING_
	if(Verbosity() >= 1) _t_seeds_cleanup->restart();
	// we still need to update the track fields for DCA and PCA
	// we can do that from the final vertex position

	shift_dx = -_vertex[0];
	shift_dy = -_vertex[1];
	shift_dz = -_vertex[2];

	// shift to precision final vertex
	shift_coordinate_system(shift_dx, shift_dy, shift_dz);

	// recompute track fits to fill dca and pca + error fields
	std::vector<SimpleTrack3D> refit_tracks;
	std::vector<double> refit_errors;
	std::vector<Eigen::Matrix<float, 5, 5> > refit_covars;

	if(Verbosity() >= 1){
	  cout<<__LINE__<< ": Event: "<< _event << ": # tracks before cleanup: "<< _tracks.size() <<endl;
	}

	_tracker->finalize(_tracks, refit_tracks);

	if(Verbosity() >= 1){
	  cout<<__LINE__<< ": Event: "<< _event << ": # tracks after cleanup: "<< _tracks.size()  <<endl;
	}

	for (unsigned int tt = 0; tt < refit_tracks.size(); ++tt) {
		refit_errors.push_back(_tracker->getKalmanStates()[tt].chi2);
		refit_covars.push_back(_tracker->getKalmanStates()[tt].C);
	}


	_tracks = refit_tracks;
	_track_errors = refit_errors;
	_track_covars = refit_covars;

	// shift back to global coordinates
	shift_coordinate_system(-shift_dx, -shift_dy, -shift_dz);
	if(Verbosity() >= 1) _t_seeds_cleanup->stop();
#endif

	// okay now we are done with the tracker
	_tracker->clear();

	//FIXME yuhw
	_clusters.clear();

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::export_output() {

	if (_all_tracks.empty())
		return Fun4AllReturnCodes::EVENT_OK;

	SvtxVertex_v1 vertex;
	vertex.set_t0(0.0);
	for (int i = 0; i < 3; ++i)
		vertex.set_position(i, _vertex[i]);
	vertex.set_chisq(0.0);
	vertex.set_ndof(0);
	vertex.set_error(0, 0, 0.0);
	vertex.set_error(0, 1, 0.0);
	vertex.set_error(0, 2, 0.0);
	vertex.set_error(1, 0, 0.0);
	vertex.set_error(1, 1, 0.0);
	vertex.set_error(1, 2, 0.0);
	vertex.set_error(2, 0, 0.0);
	vertex.set_error(2, 1, 0.0);
	vertex.set_error(2, 2, 0.0);

	// at this point we should already have an initial pt and pz guess...
	// need to translate this into the PHG4Track object...

	vector<SimpleHit3D> track_hits;
	int clusterID;

	for (unsigned int itrack = 0; itrack < _all_tracks.size(); itrack++) {
		SvtxTrack_v1 track;
		track.set_id(itrack);
		track_hits.clear();
		track_hits = _all_tracks.at(itrack).hits;

		for (unsigned int ihit = 0; ihit < track_hits.size(); ihit++) {
			if ((track_hits.at(ihit).get_id()) >= _cluster_map->size()) {
				continue;
			}
			SvtxCluster* cluster = _cluster_map->get(
					track_hits.at(ihit).get_id());
			//mark hit asu used by iteration number n
			//_hit_used_map[track_hits.at(ihit).get_id()] = _n_iteration;
			_assoc_container->SetClusterTrackAssoc(track_hits.at(ihit).get_id(), track.get_id());
			clusterID = cluster->get_id();
#ifdef _DEBUG_
			cout
			<<__LINE__
			<<": itrack: " << itrack
			<<": nhits: " << track_hits.size()
			<<": hitID: " << clusterID
			<<": layer: " << cluster->get_layer()
			<<endl;
#endif

			//TODO verify this change
			//int clusterLayer = cluster->get_layer();
			//if ((clusterLayer < (int) _nlayers_seeding) && (clusterLayer >= 0)) {
			track.insert_cluster(clusterID);
			//}
		}

		float kappa = _all_tracks.at(itrack).kappa;
		float d = _all_tracks.at(itrack).d;
		float phi = _all_tracks.at(itrack).phi;
		float dzdl = _all_tracks.at(itrack).dzdl;
		float z0 = _all_tracks.at(itrack).z0;

		//    track.set_helix_phi(phi);
		//    track.set_helix_kappa(kappa);
		//    track.set_helix_d(d);
		//    track.set_helix_z0(z0);
		//    track.set_helix_dzdl(dzdl);

		float pT = kappaToPt(kappa);

		float x_center = cos(phi) * (d + 1 / kappa); // x coordinate of circle center
		float y_center = sin(phi) * (d + 1 / kappa);  // y    "      "     " "

		// find helicity from cross product sign
		short int helicity;
		if ((track_hits[0].get_x() - x_center)
				* (track_hits[track_hits.size() - 1].get_y() - y_center)
				- (track_hits[0].get_y() - y_center)
						* (track_hits[track_hits.size() - 1].get_x() - x_center)
				> 0) {
			helicity = 1;
		} else {
			helicity = -1;
		}
		float pZ = 0;
		if (dzdl != 1) {
			pZ = pT * dzdl / sqrt(1.0 - dzdl * dzdl);
		}
		int ndf = 2 * _all_tracks.at(itrack).hits.size() - 5;
		track.set_chisq(_all_track_errors[itrack]);
		track.set_ndf(ndf);
		track.set_px(pT * cos(phi - helicity * M_PI / 2));
		track.set_py(pT * sin(phi - helicity * M_PI / 2));
		track.set_pz(pZ);

		track.set_dca2d(d);
		track.set_dca2d_error(sqrt(_all_track_covars[itrack](1, 1)));

		if (_magField > 0) {
			track.set_charge(helicity);
		} else {
			track.set_charge(-1.0 * helicity);
		}

		Eigen::Matrix<float, 6, 6> euclidean_cov =
				Eigen::Matrix<float, 6, 6>::Zero(6, 6);
		convertHelixCovarianceToEuclideanCovariance(_magField, phi, d, kappa,
				z0, dzdl, _all_track_covars[itrack], euclidean_cov);

		for (unsigned int row = 0; row < 6; ++row) {
			for (unsigned int col = 0; col < 6; ++col) {
				track.set_error(row, col, euclidean_cov(row, col));
			}
		}

		track.set_x(vertex.get_x() + d * cos(phi));
		track.set_y(vertex.get_y() + d * sin(phi));
		track.set_z(vertex.get_z() + z0);

		_track_map->insert(&track);
		vertex.insert_track(track.get_id());

		if (Verbosity() > 5) {
			cout << "track " << itrack << " quality = " << track.get_quality()
					<< endl;
			cout << "px = " << track.get_px() << " py = " << track.get_py()
					<< " pz = " << track.get_pz() << endl;
		}
	}  // track loop

	SvtxVertex *vtxptr = _vertex_map->insert(&vertex);
	if (Verbosity() > 5)
		vtxptr->identify();

	if (Verbosity() > 0) {
		cout << "PHHoughSeeding::process_event -- leaving process_event"
				<< endl;
	}

	// we are done with these now...
	_clusters.clear();
	_all_tracks.clear();
	_all_track_errors.clear();
	_all_track_covars.clear();
	_vertex.clear();
	_vertex.assign(3, 0.0);

	return Fun4AllReturnCodes::EVENT_OK;
}

float PHHoughSeeding::kappaToPt(float kappa) {
	return _pt_rescale * _magField / 333.6 / kappa;
}

float PHHoughSeeding::ptToKappa(float pt) {
	return _pt_rescale * _magField / 333.6 / pt;
}

void PHHoughSeeding::convertHelixCovarianceToEuclideanCovariance(float B,
		float phi, float d, float kappa, float z0, float dzdl,
		Eigen::Matrix<float, 5, 5> const& input,
		Eigen::Matrix<float, 6, 6>& output) {

	Eigen::Matrix<float, 6, 5> J = Eigen::Matrix<float, 6, 5>::Zero(6, 5);

	// phi,d,nu,z0,dzdl
	// -->
	// x,y,z,px,py,pz

	float nu = sqrt(kappa);
	float dk_dnu = 2 * nu;

	float cosphi = cos(phi);
	float sinphi = sin(phi);

	J(0, 0) = -sinphi * d;
	J(0, 1) = cosphi;
	J(1, 0) = d * cosphi;
	J(1, 1) = sinphi;
	J(2, 3) = 1;

	float pt = 0.003 * B / kappa;
	float dpt_dk = -0.003 * B / (kappa * kappa);

	J(3, 0) = -cosphi * pt;
	J(3, 2) = -sinphi * dpt_dk * dk_dnu;
	J(4, 0) = -sinphi * pt;
	J(4, 2) = cosphi * dpt_dk * dk_dnu;

	float alpha = 1. / (1. - dzdl * dzdl);
	float alpha_half = pow(alpha, 0.5);
	float alpha_3_half = alpha * alpha_half;

	J(5, 2) = dpt_dk * dzdl * alpha_half * dk_dnu;
	J(5, 4) = pt * (alpha_half + dzdl * dzdl * alpha_3_half) * dk_dnu;

	output = J * input * (J.transpose());
}

void PHHoughSeeding::shift_coordinate_system(double dx, double dy,
		double dz) {

	for (unsigned int ht = 0; ht < _clusters.size(); ++ht) {
		_clusters[ht].set_x(_clusters[ht].get_x() + dx);
		_clusters[ht].set_y(_clusters[ht].get_y() + dy);
		_clusters[ht].set_z(_clusters[ht].get_z() + dz);
	}

	for (unsigned int tt = 0; tt < _tracks.size(); ++tt) {
		for (unsigned int hh = 0; hh < _tracks[tt].hits.size(); ++hh) {
			_tracks[tt].hits[hh].set_x(_tracks[tt].hits[hh].get_x() + dx);
			_tracks[tt].hits[hh].set_y(_tracks[tt].hits[hh].get_y() + dy);
			_tracks[tt].hits[hh].set_z(_tracks[tt].hits[hh].get_z() + dz);
		}
	}

	_vertex[0] += dx;
	_vertex[1] += dy;
	_vertex[2] += dz;

	return;
}

int PHHoughSeeding::CleanupSeedsByHitPattern() {

        std::vector<SimpleTrack3D> _tracks_cleanup;
        _tracks_cleanup.clear();

	if(Verbosity() >= 1){
	  cout<<__LINE__<< ": Event: "<< _event << ": # tracks before cleanup: "<< _tracks.size() <<endl;
	}

	/*
	  std::vector<double> _track_errors_cleanup;
	  _track_errors_cleanup.clear();	  
	  std::vector<Eigen::Matrix<float, 5, 5> > _track_covars_cleanup;
	  _track_covars_cleanup.clear();
	  
	  std::vector<HelixKalmanState> _kalman_states_cleanup;
	  _kalman_states_cleanup.clear();
	*/
       
	typedef std::tuple<int, int, int, int> KeyType;
	typedef std::multimap< KeyType, unsigned int > MapKeyTrkID;
	
	std::set<KeyType> keys;
	std::vector<bool> v_track_used;
	MapKeyTrkID m_key_itrack;


	typedef std::set<unsigned int> TrackList;

	std::set<unsigned int> OutputList;
	std::multimap<int, unsigned int > hitIdTrackList;

	unsigned int max_hit_id = 0;
	//For each hit make list of all associated tracks

	std::vector<bool> good_track;
	//	printf("build hit track map\n");
	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
	  good_track.push_back(true);
	  SimpleTrack3D track = _tracks[itrack];
	  for( SimpleHit3D hit : track.hits) {
	    hitIdTrackList.insert(std::make_pair(hit.get_id(),itrack));
	    if(hit.get_id()>max_hit_id) max_hit_id = hit.get_id();
	  }
	}
	//	printf("build track duplicate map\n");
	//Check Tracks for duplicates by looking for hits shared
	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
	  if(good_track[itrack]==false) continue;//already checked this one
	  if(OutputList.count(itrack)>0)continue;//already got this one
	  
	  SimpleTrack3D track = _tracks[itrack];

	  int trackid_max_nhit = itrack;
	  unsigned int max_nhit = track.hits.size();
	  int onhit = track.hits.size();

	  TrackList tList;
	  for( SimpleHit3D hit : track.hits) {
	    int nmatch = hitIdTrackList.count(hit.get_id());
	    if(nmatch>1){
	      multimap<int, unsigned int >::iterator it = hitIdTrackList.find(hit.get_id());
	      //Loop over track matches and add them to the list, select longest in the process
	      for(; it != hitIdTrackList.end();++it) {
		unsigned int match_trackid = (*it).second;
		if(match_trackid == itrack) continue;//original track
		if(good_track[match_trackid]==false)continue;
		tList.insert(match_trackid);
		SimpleTrack3D mtrack = _tracks[match_trackid];
	      }
	    }
	  }
	  //	  int tlsize = tList.size();

	  //	  cout << "remove bad matches " << tList.size() << "itrk: " << itrack << endl;
	  //loop over matches and remove matches with too few shared hits  
	  TrackList mergeList;
	  for( unsigned int match : tList) {
	    //	    cout << "processing " << match << " of " << tList.size() << " itrk " << itrack << endl;
	    if(match==itrack)continue;
	    if(good_track[match]==false)continue;

	    SimpleTrack3D mtrack = _tracks[match]; //matched track
	    int mnhit = mtrack.hits.size();
 	    std::set<unsigned int> HitList;
	    //put hits from both tracks in a set
	    for( SimpleHit3D hit : track.hits) HitList.insert(hit.get_id());
	    for( SimpleHit3D hit : mtrack.hits) HitList.insert(hit.get_id());
	    //set stores only unique hits, tracks overlap if:
	    int sumnhit = HitList.size();
	    if(sumnhit<(onhit+mnhit-3)){// more than 3 overlaps 
	      //not enough overlap, drop track from list
	      //tList.erase(match);
	      //good_track[match] = false;
	      if(sumnhit != onhit){//no subset
		mergeList.insert(match);
	      }
	    }

	  }

	  tList.clear();
	  //	  cout << "flag bad matches done " << mergeList.size() << " itrk " << itrack << endl;
	  //loop over matches and flag all tracks bad except the longest 
	  std::set<unsigned int> MergedHitList;
	  if(mergeList.size()==0){
	    for( SimpleHit3D hit : track.hits) MergedHitList.insert(hit.get_id());
	  }
	  //	  cout << "merge good matches itrk " << itrack << " #" << mergeList.size() << endl;
	  for( unsigned int match : mergeList) {
	    if(match==itrack)continue;
	    if(good_track[match]==false)continue;
	    //	    cout << "  adding " << match << endl;
	    //check number of shared hits
	    //get tracks

	    SimpleTrack3D mtrack = _tracks[match]; //matched track
	    if(mtrack.hits.size()>max_nhit){
	      max_nhit = mtrack.hits.size();
	      trackid_max_nhit = match;
	      good_track[itrack] = false;
	    }else{
	      good_track[match] = false;
	    }
	    for( SimpleHit3D hit : track.hits) MergedHitList.insert(hit.get_id());
	    for( SimpleHit3D hit : mtrack.hits) MergedHitList.insert(hit.get_id());
	  }

	  //	  int ntracks = _tracks.size();
	  //int outtracks = OutputList.size();
	  //	  printf("CLEANUP: itrack: %5d(%d) => %5d matches max %d(%d) tracks kept: %d\n",
	  //	 itrack, ntracks,tlsize, max_nhit, trackid_max_nhit, outtracks);

	  //	  printf("keep track %d\n",trackid_max_nhit);
	  //add merged hit list to merged track 
	  if(OutputList.count(trackid_max_nhit)==0){
	    _tracks_cleanup.push_back(_tracks[trackid_max_nhit]);
	
	    /*  _kalman_states_cleanup.push_back((_tracker->getKalmanStates())[trackid_max_nhit]);
		_track_covars_cleanup.push_back((_tracker->getKalmanStates())[trackid_max_nhit].C);
		_track_errors_cleanup.push_back(_tracker->getKalmanStates()[trackid_max_nhit].chi2);
	    */
	  }
	  OutputList.insert(trackid_max_nhit);
	  
	  _tracks_cleanup.back().hits.clear();
	  
	  for(unsigned int hitID : MergedHitList) {
	    SimpleHit3D hit;
	    hit.set_id(hitID);
	    _tracks_cleanup.back().hits.push_back(hit);
	  }
	  
				    
	}

	_tracks.clear();
	_tracks = _tracks_cleanup;

	/*	_track_errors.clear();
	  _track_errors = _track_errors_cleanup;
	  _track_covars.clear();
	  _track_covars = _track_covars_cleanup;
	  _tracker->getKalmanStates().clear();
	  for(auto &kstate :  _kalman_states_cleanup){
	  _tracker->getKalmanStates().push_back(kstate);
	  }
	*/
		
	if(Verbosity() >= 1){
	  cout<<__LINE__<< ": Event: "<< _event <<endl;
	  cout << ": # tracks after cleanup: "<< _tracks.size() << " ol:" <<OutputList.size()  <<endl;
	}

	
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHHoughSeeding::CleanupTracksByHitPattern() {

        std::vector<SimpleTrack3D> _tracks_cleanup;
        _tracks_cleanup.clear();

	//	if(Verbosity() >= 1)
	{
	    cout<<__LINE__<< ": Event: "<< _event << ": # tracks before cleanup: "<< _tracks.size() <<endl;
	  }

	
	std::vector<double> _track_errors_cleanup;
	_track_errors_cleanup.clear();	  
	std::vector<Eigen::Matrix<float, 5, 5> > _track_covars_cleanup;
	_track_covars_cleanup.clear();
	
	std::vector<HelixKalmanState> _kalman_states_cleanup;
	_kalman_states_cleanup.clear();
	       
	typedef std::tuple<int, int, int, int> KeyType;
	typedef std::multimap< KeyType, unsigned int > MapKeyTrkID;
	
	std::set<KeyType> keys;
	std::vector<bool> v_track_used;
	MapKeyTrkID m_key_itrack;


	typedef std::set<unsigned int> TrackList;

	std::set<unsigned int> OutputList;
	std::multimap<int, unsigned int > hitIdTrackList;

	unsigned int max_hit_id = 0;
	//For each hit make list of all associated tracks

	std::vector<bool> good_track;
	//	printf("build hit track map\n");
	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
	  good_track.push_back(true);
	  SimpleTrack3D track = _tracks[itrack];
	  for( SimpleHit3D hit : track.hits) {
	    hitIdTrackList.insert(std::make_pair(hit.get_id(),itrack));
	    if(hit.get_id()>max_hit_id) max_hit_id = hit.get_id();
	  }
	}
	//	printf("build track duplicate map\n");
	//Check Tracks for duplicates by looking for hits shared
	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
	  if(good_track[itrack]==false) continue;//already checked this one
	  if(OutputList.count(itrack)>0)continue;//already got this one
	  
	  SimpleTrack3D track = _tracks[itrack];

	  int trackid_max_nhit = itrack;
	  unsigned int max_nhit = track.hits.size();
	  int onhit = track.hits.size();

	  TrackList tList;
	  for( SimpleHit3D hit : track.hits) {
	    int nmatch = hitIdTrackList.count(hit.get_id());
	    if(nmatch>1){
	      multimap<int, unsigned int >::iterator it = hitIdTrackList.find(hit.get_id());
	      //Loop over track matches and add them to the list, select longest in the process
	      for(; it != hitIdTrackList.end();++it) {
		unsigned int match_trackid = (*it).second;
		if(match_trackid == itrack) continue;//original track
		if(good_track[match_trackid]==false)continue;
		tList.insert(match_trackid);
		SimpleTrack3D mtrack = _tracks[match_trackid];
	      }
	    }
	  }
	  //	  int tlsize = tList.size();

	  //	  cout << "remove bad matches " << tList.size() << "itrk: " << itrack << endl;
	  //loop over matches and remove matches with too few shared hits  
	  TrackList mergeList;
	  for( unsigned int match : tList) {
	    //	    cout << "processing " << match << " of " << tList.size() << " itrk " << itrack << endl;
	    if(match==itrack)continue;
	    if(good_track[match]==false)continue;

	    SimpleTrack3D mtrack = _tracks[match]; //matched track
	    int mnhit = mtrack.hits.size();
 	    std::set<unsigned int> HitList;
	    //put hits from both tracks in a set
	    for( SimpleHit3D hit : track.hits) HitList.insert(hit.get_id());
	    for( SimpleHit3D hit : mtrack.hits) HitList.insert(hit.get_id());
	    //set stores only unique hits, tracks overlap if:
	    int sumnhit = HitList.size();
	    if(sumnhit<(onhit+mnhit-3)){// more than 3 overlaps 
	      //not enough overlap, drop track from list
	      //tList.erase(match);
	      //good_track[match] = false;
	      if(sumnhit != onhit){//no subset
		mergeList.insert(match);
	      }
	    }

	  }

	  tList.clear();
	  //	  cout << "flag bad matches done " << mergeList.size() << " itrk " << itrack << endl;
	  //loop over matches and flag all tracks bad except the longest 
	  std::set<unsigned int> MergedHitList;
	  if(mergeList.size()==0){
	    for( SimpleHit3D hit : track.hits) MergedHitList.insert(hit.get_id());
	  }
	  //	  cout << "merge good matches itrk " << itrack << " #" << mergeList.size() << endl;
	  for( unsigned int match : mergeList) {
	    if(match==itrack)continue;
	    if(good_track[match]==false)continue;
	    //	    cout << "  adding " << match << endl;
	    //check number of shared hits
	    //get tracks

	    SimpleTrack3D mtrack = _tracks[match]; //matched track
	    if(mtrack.hits.size()>max_nhit){
	      max_nhit = mtrack.hits.size();
	      trackid_max_nhit = match;
	      good_track[itrack] = false;
	    }else{
	      good_track[match] = false;
	    }
	    for( SimpleHit3D hit : track.hits) MergedHitList.insert(hit.get_id());
	    for( SimpleHit3D hit : mtrack.hits) MergedHitList.insert(hit.get_id());
	  }

	  //	  int ntracks = _tracks.size();
	  //int outtracks = OutputList.size();
	  //	  printf("CLEANUP: itrack: %5d(%d) => %5d matches max %d(%d) tracks kept: %d\n",
	  //	 itrack, ntracks,tlsize, max_nhit, trackid_max_nhit, outtracks);

	  //	  printf("keep track %d\n",trackid_max_nhit);
	  //add merged hit list to merged track 
	  if(OutputList.count(trackid_max_nhit)==0){
	    _tracks_cleanup.push_back(_tracks[trackid_max_nhit]);
	
	    _kalman_states_cleanup.push_back((_tracker->getKalmanStates())[trackid_max_nhit]);
	    _track_covars_cleanup.push_back((_tracker->getKalmanStates())[trackid_max_nhit].C);
	    _track_errors_cleanup.push_back(_tracker->getKalmanStates()[trackid_max_nhit].chi2);

	  }
	  OutputList.insert(trackid_max_nhit);
	  
	  _tracks_cleanup.back().hits.clear();
	  
	  for(unsigned int hitID : MergedHitList) {
	    SimpleHit3D hit;
	    hit.set_id(hitID);
	    _tracks_cleanup.back().hits.push_back(hit);
	  }
	  
				    
	}

	_tracks.clear();
	_tracks = _tracks_cleanup;

	_track_errors.clear();
	_track_errors = _track_errors_cleanup;
	_track_covars.clear();
	_track_covars = _track_covars_cleanup;
	_tracker->getKalmanStates().clear();
	for(auto &kstate :  _kalman_states_cleanup){
	  _tracker->getKalmanStates().push_back(kstate);
	}
	
		
       	if(Verbosity() >= 1)
	  {
	    cout<<__LINE__<< ": Event: "<< _event <<endl;
	    cout << ": # tracks after cleanup: "<< _tracks.size() << " ol:" <<OutputList.size()  <<endl;
	  }
	
	return Fun4AllReturnCodes::EVENT_OK;
}


int PHHoughSeeding::CleanupSeeds() {

	std::vector<SimpleTrack3D> _tracks_cleanup;
	_tracks_cleanup.clear();

	typedef std::tuple<int, int, int, int> KeyType;
	typedef std::multimap< KeyType, unsigned int > MapKeyTrkID;

	std::set<KeyType> keys;
	std::vector<bool> v_track_used;
	MapKeyTrkID m_key_itrack;

#ifdef _DEBUG_
	cout<<__LINE__<<": CleanupSeeds: Event: "<<_event<<endl;
#endif

	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
		SimpleTrack3D track = _tracks[itrack];

		int id = track.d / _max_merging_dr;
		int iz = track.z0 / _max_merging_dz;
		int iphi = track.phi / _max_merging_dphi;
		int idzdl = track.dzdl / _max_merging_deta;

#ifdef _DEBUG_
//		printf("itrack: %d: \t%e, \t%e, \t%e, \t%e, %5d, %5d, %5d, %5d \n",
//				itrack,
//				track.d, track.z0, track.phi, track.dzdl,
//				id, iz, iphi, idzdl
//		);
#endif

		KeyType key = std::make_tuple(id, iz, iphi, idzdl);

		keys.insert(key);
		m_key_itrack.insert(
				std::make_pair(key,
						itrack));

		v_track_used.push_back(false);
	}

#ifdef _DEBUG_
	for(auto it = m_key_itrack.begin();
						it != m_key_itrack.end();
						++it) {

		KeyType key = it->first;
		unsigned int itrack = it->second;

		int id = std::get < 0 > (key);
		int iz = std::get < 1 > (key);
		int iphi = std::get < 2 > (key);
		int idzdl = std::get < 3 > (key);

		SimpleTrack3D track = _tracks[itrack];

		cout << __LINE__ << endl;
		printf("itrack: %5d => {%5d, %5d, %5d, %5d} \n",
				itrack,
				id, iz, iphi, idzdl);
	}
#endif

	for(KeyType key : keys) {

		unsigned int itrack = m_key_itrack.equal_range(key).first->second;

#ifdef _DEBUG_
		cout<<"---------------------------------------------------\n";
		cout<<__LINE__<< ": processing: "<<itrack<<endl;
		cout<<"---------------------------------------------------\n";
#endif

		if (v_track_used[itrack] == true)
			continue;
#ifdef _DEBUG_
		cout<<__LINE__<<":    itrack: "<< itrack <<": {";
#endif
		std::set<unsigned int> hitIDs;
		for( SimpleHit3D hit : _tracks[itrack].hits) {
			hitIDs.insert(hit.get_id());
#ifdef _DEBUG_
			cout<<hit.get_id() <<", ";
#endif
		}
#ifdef _DEBUG_
		cout<<"}"<<endl;
#endif

		//! find tracks winthin neighbor bins
		std::vector<unsigned int> v_related_tracks;
		for (int id = std::get < 0 > (key) - 1; id <= std::get < 0 > (key) + 1;
				++id) {
			for (int iz = std::get < 1 > (key) - 1;
					iz <= std::get < 1 > (key) + 1; ++iz) {
				for (int iphi = std::get < 2 > (key) - 1;
						iphi <= std::get < 2 > (key) + 1; ++iphi) {
					for (int idzdl = std::get < 3 > (key) - 1;
							idzdl <= std::get < 3 > (key) + 1; ++idzdl) {
						KeyType key_temp = std::make_tuple(id, iz, iphi, idzdl);

						if (m_key_itrack.find(key_temp) != m_key_itrack.end()) {
							for (auto it =
									m_key_itrack.equal_range(key_temp).first;
									it
											!= m_key_itrack.equal_range(
													key_temp).second; ++it) {

								if(it->second == itrack)
									continue;

								unsigned int share_hits = 0;
								for(SimpleHit3D hit : _tracks[it->second].hits) {
									unsigned int hitID = hit.get_id();
									if(std::find(
											hitIDs.begin(),
											hitIDs.end(),
											hitID) != hitIDs.end()) {

										++share_hits;
										if(share_hits > _max_share_hits) {
											v_related_tracks.push_back(it->second);
#ifdef _DEBUG_
											cout<<__LINE__<<": rel track: "<<it->second <<": {";
											for(SimpleHit3D hit : _tracks[it->second].hits) {
												cout<< hit.get_id() <<", ";
											}
											cout<<"}"<<endl;
#endif
											break;
										}
									}
								} //loop to find common hits

							}
//#ifdef _DEBUG_
//							cout << __LINE__ << ": ";
//							printf("{%5d, %5d, %5d, %5d} => {%d, %d} \n", id,
//									iz, iphi, idzdl,
//									m_key_itrack.equal_range(key_temp).first->second,
//									m_key_itrack.equal_range(key_temp).second->second);
//#endif
						}
					}
				}
			}
		}

		if(v_related_tracks.size() == 0) {
			_tracks_cleanup.push_back(_tracks[itrack]);
		} else {

			_tracks_cleanup.push_back(_tracks[itrack]);

			_tracks_cleanup.back().hits.clear();

#ifdef _DEBUG_
			int n_merge_track = 1;
			cout<<__LINE__<<": nclusters before merge: "<< hitIDs.size() <<endl;
#endif

			//! Add hits from other related tracks
			//std::set<unsigned int> hitIDs;
			for(unsigned int irel : v_related_tracks) {

				if(v_track_used[irel] == true) continue;

				//! hits from itrack already registered
				//if(irel == itrack) continue;

#ifdef _DEBUG_
				++n_merge_track;
#endif

#ifdef _MEARGE_SEED_CLUSTER_
				SimpleTrack3D track = _tracks[irel];
				for(SimpleHit3D hit : track.hits) {
					hitIDs.insert(hit.get_id());
				}
#endif
				v_track_used[irel] = true;
			}

#ifdef _DEBUG_
			cout<<__LINE__<<": # tracks merged: "<< n_merge_track <<endl;
			cout<<"{ ";
#endif
			for(unsigned int hitID : hitIDs) {
				SimpleHit3D hit;
				hit.set_id(hitID);
#ifdef _DEBUG_
				cout<<hitID <<", ";
#endif
				_tracks_cleanup.back().hits.push_back(hit);
			}
#ifdef _DEBUG_
			cout<<"}"<<endl;
			cout<<__LINE__<<": nclusters after merge:  "<< hitIDs.size() <<endl;
			cout<<__LINE__<<": nclusters after merge:  "<< _tracks_cleanup.back().hits.size() <<endl;
#endif
		}

		v_track_used[itrack] = true;
	}

#ifdef _DEBUG_
	cout<<__LINE__<< ": Event: "<< _event <<endl;
	cout << ": # tracks before cleanup: "<< _tracks.size() <<endl;
	cout << ": # tracks after  cleanup: "<< _tracks_cleanup.size() <<endl;
#endif
	_tracks.clear();
	_tracks = _tracks_cleanup;

	return Fun4AllReturnCodes::EVENT_OK;
}



/*!
 *  \file PHG4KalmanPatRec.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#include "PHG4KalmanPatRec.h"

// g4hough includes
#include "SvtxVertexMap.h"
#include "SvtxVertexMap_v1.h"
#include "SvtxVertex.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxTrackState.h"
#include "SvtxClusterMap.h"
#include "SvtxCluster.h"
#include "SvtxHit_v1.h"
#include "SvtxHitMap.h"

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeomSiLadders.h>

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

PHG4KalmanPatRec::PHG4KalmanPatRec(
		const string& name,
		unsigned int nlayers_maps,
		unsigned int nlayers_intt,
		unsigned int nlayers_tpc,
		unsigned int nlayers_seeding,
		unsigned int min_nlayers_seeding
		)
    : SubsysReco(name),
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
      _g4clusters(NULL),
      _g4tracks(NULL),
      _g4vertexes(NULL),
      _svtxhitsmap(nullptr),
      _hit_used_map(NULL),
      _cells_svtx(nullptr),
      _cells_intt(nullptr),
      _cells_maps(nullptr),
      _geom_container_intt(nullptr),
      _geom_container_maps(nullptr),
      _n_iteration(0),
      _n_max_iterations(2),
      _seeding_only_mode(false),
      _analyzing_mode(false),
      _analyzing_file(NULL),
      _analyzing_ntuple(NULL),
      _max_merging_dphi(0.1),
      _max_merging_deta(0.1),
      _max_merging_dr(0.1),
      _max_merging_dz(0.1),
      _max_share_hits(3),
      _fitter(NULL),
      //      _track_fitting_alg_name("DafRef"),
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
  
  _user_material.clear();
	for(unsigned int i=0;i<_nlayers_maps;++i)
		_user_material[i] = 0.003;
	for(unsigned int i=_nlayers_maps;i<_nlayers_maps+_nlayers_intt;++i)
		_user_material[i] = 0.008;

//	unsigned int maps_layers[] = {0, 1, 2};
//	this->set_maps_layers(maps_layers, 3);
//
//	unsigned int intt_layers[] = {3, 4, 5, 6};
//	this->set_intt_layers(intt_layers, 4);

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
}

int PHG4KalmanPatRec::Init(PHCompositeNode* topNode) {
	if(_analyzing_mode){
	  cout << "Ana Mode, creating ntuples! " << endl;
	  _analyzing_file = new TFile("./PatRecAnalysis.root","RECREATE");
	  //	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","spt:seta:sphi:pt:eta:phi:layer:ncand:nmeas");
	  _analyzing_ntuple = new TNtuple("ana_nt","ana_nt","pt:kappa:d:phi:dzdl:z0:nhit:ml:rec:dt");
	  cout << "Done" << endl;
	  
	}
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4KalmanPatRec::InitRun(PHCompositeNode* topNode) {

	int code = Fun4AllReturnCodes::ABORTRUN;

	code = CreateNodes(topNode);
	if(code != Fun4AllReturnCodes::EVENT_OK)
		return code;

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
	
	code = InitializeGeometry(topNode);
	if(code != Fun4AllReturnCodes::EVENT_OK)
	  return code;
	code = InitializePHGenFit(topNode);
	if(code != Fun4AllReturnCodes::EVENT_OK)
		return code;

	/*!
	 * Initilize parameters
	 */
	for(int layer = 0; layer < _nlayers_all; ++layer) {
		_search_wins_phi.insert(std::make_pair(layer, _search_win_phi));
		_search_wins_theta.insert(std::make_pair(layer, _search_win_theta));
		_max_incr_chi2s.insert(std::make_pair(layer, _max_incr_chi2));
	}

	// nightly build 2017-05-04
	//	_search_wins_phi[8]  = 50.;
	//	_search_wins_phi[9]  = 45.;
	//	_search_wins_phi[10] = 40.;
//	_search_wins_phi[11] = 30.;
//	_search_wins_phi[12] = 30.;
//	_search_wins_phi[13] = 30.;
//	_search_wins_phi[14] = 30.;
//	_search_wins_phi[15] = 30.;
//	_search_wins_phi[16] = 30.;
//	_search_wins_phi[17] = 30.;
//	_search_wins_phi[18] = 30.;
//	_search_wins_phi[19] = 30.;
//	_search_wins_phi[20] = 30.;
//
//	_max_incr_chi2s[8]  = _max_incr_chi2s[8] < 1000. ? 1000 : _max_incr_chi2s[8];
//	_max_incr_chi2s[9]  = _max_incr_chi2s[9] < 500.  ? 500  : _max_incr_chi2s[9];
//	_max_incr_chi2s[10] = _max_incr_chi2s[10]< 500.  ? 500  : _max_incr_chi2s[10];
//	_max_incr_chi2s[11] = _max_incr_chi2s[11]< 200.  ? 200  : _max_incr_chi2s[11];
//	_max_incr_chi2s[12] = _max_incr_chi2s[12]< 200.  ? 200  : _max_incr_chi2s[12];
//	_max_incr_chi2s[13] = _max_incr_chi2s[13]< 100.  ? 100  : _max_incr_chi2s[13];
//	_max_incr_chi2s[14] = _max_incr_chi2s[14]< 100.  ? 100  : _max_incr_chi2s[14];
//	_max_incr_chi2s[15] = _max_incr_chi2s[15]< 100.  ? 100  : _max_incr_chi2s[15];
//	_max_incr_chi2s[16] = _max_incr_chi2s[16]< 100.  ? 100  : _max_incr_chi2s[16];
//	_max_incr_chi2s[17] = _max_incr_chi2s[17]< 100.  ? 100  : _max_incr_chi2s[17];
//	_max_incr_chi2s[18] = _max_incr_chi2s[18]< 50.   ? 50   : _max_incr_chi2s[18];
//	_max_incr_chi2s[19] = _max_incr_chi2s[19]< 50.   ? 50   : _max_incr_chi2s[19];
//	_max_incr_chi2s[20] = _max_incr_chi2s[20]< 50.   ? 50   : _max_incr_chi2s[20];

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
				<< "====================== PHG4KalmanPatRec::InitRun() ======================"
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
	
	return code;
}

int PHG4KalmanPatRec::process_event(PHCompositeNode *topNode) {

  if (Verbosity() > 0){
	  cout << "PHG4KalmanPatRec::process_event -- entered" << endl;
	  cout << "nMapsLayers = " << _nlayers_maps << endl;
	  cout << "nInttLayers = " << _nlayers_intt << endl;
	  cout << "nTPCLayers = " << _nlayers_tpc << endl;
  }
	// start fresh
	int code;
	_n_iteration = 0;
	if(_n_max_iterations<1)_n_max_iterations = 1;
	if(_n_max_iterations>3)_n_max_iterations = 3;
	_clusters.clear();
	_all_tracks.clear();
	_all_track_errors.clear();
	_all_track_covars.clear();

	_vertex.clear();
	_vertex.assign(3, 0.0);
	//-----------------------------------
	// Get Objects off of the Node Tree
	//-----------------------------------

	GetNodes(topNode);// Allocate Cluster Use Map allocated in here
	
	for(_n_iteration = 1;_n_iteration<=_n_max_iterations;_n_iteration++){
	  _tracks.clear();
	  _track_errors.clear();
	  _track_covars.clear();
	  
	  if(_n_iteration==1){    
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
	    code = InitializeGeometry(topNode);
	    if(Verbosity() >= 1) _t_seed_init1->restart();
	    if(code != Fun4AllReturnCodes::EVENT_OK)
	      return code;
	  }
	  
	  if(_n_iteration==2){
	    int min_layers    = 6;
	    int nlayers_seeds = 12;
	    int seeding_layers[] = {(int)(_nlayers_maps+_nlayers_intt),
				    (int)(_nlayers_maps+_nlayers_intt+1),
				    (int)(_nlayers_maps+_nlayers_intt+2),
				    
				    (int)(_nlayers_maps+_nlayers_intt+9),
				    (int)(_nlayers_maps+_nlayers_intt+10),
				    (int)(_nlayers_maps+_nlayers_intt+11),

				    (int)(_nlayers_maps+_nlayers_intt+20),
				    (int)(_nlayers_maps+_nlayers_intt+21),

				    (int)(_nlayers_maps+_nlayers_intt+31),
				    (int)(_nlayers_maps+_nlayers_intt+32),

				    (int)(_nlayers_maps+_nlayers_intt+38),

				    (int)(_nlayers_maps+_nlayers_intt+45)  // avoid the outer TPC layer, it is inefficient
				    //7,13,19,25,31,37,46
				    //7,8,13,14,19,20,26,27,34,35,40,46
	    };
	    
	    set_seeding_layer(seeding_layers, nlayers_seeds);
	    set_min_nlayers_seeding(min_layers);
	    _min_combo_hits = min_layers;
	    _max_combo_hits = nlayers_seeds;
	    code = InitializeGeometry(topNode);
	    if(Verbosity() >= 1) _t_seed_init2->restart();

	    if(code != Fun4AllReturnCodes::EVENT_OK)
	      return code;
	  }
	  if(_n_iteration==3){
	    int min_layers    = 4;
	    int nlayers_seeds = 12;
	    
	    int seeding_layers[] = {(int)(_nlayers_maps+_nlayers_intt),
				    (int)(_nlayers_maps+_nlayers_intt+1),
				    (int)(_nlayers_maps+_nlayers_intt+8),
				    (int)(_nlayers_maps+_nlayers_intt+9),
				    (int)(_nlayers_maps+_nlayers_intt+16),
				    (int)(_nlayers_maps+_nlayers_intt+17),
				    (int)(_nlayers_maps+_nlayers_intt+24),
				    (int)(_nlayers_maps+_nlayers_intt+25),
				    (int)(_nlayers_maps+_nlayers_intt+32),
				    (int)(_nlayers_maps+_nlayers_intt+33),
				    (int)(_nlayers_maps+_nlayers_intt+40),
				    (int)(_nlayers_maps+_nlayers_intt+45)  // avoid the outer TPC layer, it is inefficient
				    //7,13,19,25,31,37,46
				    //7,8,13,14,19,20,26,27,34,35,40,46
	    };
	    
	    set_seeding_layer(seeding_layers, nlayers_seeds);
	    set_min_nlayers_seeding(min_layers);
	    _min_combo_hits = min_layers;
	    _max_combo_hits = nlayers_seeds;

	    code = InitializeGeometry(topNode);
	    if(Verbosity() >= 1) _t_seed_init3->restart();
	    if(code != Fun4AllReturnCodes::EVENT_OK)
	      return code;
	  }

	  if(Verbosity() >= 1)
	    cout << "Iteration number " << _n_iteration << endl; 
	  _min_nlayers_seeding--;
	  if(Verbosity() >= 1) _t_seeding->restart();
	  
	  //-----------------------------------
	  // Translate into Helix_Hough objects
	  //-----------------------------------
	  
	  code = translate_input();//Check if cluster is already used in here
	  if (code != Fun4AllReturnCodes::EVENT_OK) return code;
	  
	  //-----------------------------------
	  // Guess a vertex position
	  //-----------------------------------
	  
	  //	code = fast_vertex_guessing();
	  //	if (code != Fun4AllReturnCodes::EVENT_OK) return code;
	  // here expect vertex to be better than +/-2.0 cm
	  
	  //-----------------------------------
	  // Find an initial vertex with tracks
	  //-----------------------------------
	  
	  //	code = initial_vertex_finding();
	  //	if (code != Fun4AllReturnCodes::EVENT_OK) return code;
	  if(_n_iteration ==1){
	    code = vertexing(topNode);
	    if (code != Fun4AllReturnCodes::EVENT_OK) return code;
	    // here expect vertex to be better than +/- 500 um
	  }
	  //-----------------------------------
	  // Seeding
	  //-----------------------------------
	  //TODO simplify this function
	  code = full_track_seeding();
	  if (code != Fun4AllReturnCodes::EVENT_OK)
	    return code;
	  
	  if(Verbosity() >= 1) _t_seeding->stop();
	  if(Verbosity() >= 1&&_n_iteration==3) _t_seed_init3->stop();
	  if(Verbosity() >= 1&&_n_iteration==2) _t_seed_init2->stop();
	  if(Verbosity() >= 1&&_n_iteration==1) _t_seed_init1->stop();

	  if(Verbosity() >= 1) _t_kalman_pat_rec->restart();
	  
	  //-----------------------------------
	  // Kalman cluster association
	  //-----------------------------------
	  if (!_seeding_only_mode) {
	    code = FullTrackFitting(topNode);
	    if (code != Fun4AllReturnCodes::EVENT_OK)
	      return code;
	  }
	  if(Verbosity() >= 1) _t_kalman_pat_rec->stop();
	  
	  //-----------------------------------
	  // Translate back into SVTX objects
	  //-----------------------------------
	  
	  //	  add_tracks();
	  
	  if(Verbosity() > 1) print_timers();
	  
	  
	}
	//	CleanupTracksByHitPattern();
       
	if(!_seeding_only_mode)
	  code = ExportOutput();
	else
	  code = export_output();
	if (code != Fun4AllReturnCodes::EVENT_OK)
	  return code;
	++_event;
	
	return Fun4AllReturnCodes::EVENT_OK;
}


void PHG4KalmanPatRec::print_timers() {
  
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

int PHG4KalmanPatRec::End(PHCompositeNode *topNode) {

	if (_do_evt_display)
		_fitter->displayEvent();

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

void PHG4KalmanPatRec::projectToRadius(const SvtxTrack* track, double B,
		double radius, std::vector<double>& intersection) {
	intersection.clear();
	intersection.assign(3, NAN);

	// start from the inner most state vector
	const SvtxTrackState* state = track->get_state(0.0);
	projectToRadius(state, track->get_charge(), B, radius, intersection);

	// iterate once to see if there is a state vector closer to the intersection
	if (track->size_states() == 1)
		return;

	const SvtxTrackState* closest = NULL;
	float min_dist = FLT_MAX;
	for (SvtxTrack::ConstStateIter iter = track->begin_states();
			iter != track->end_states(); ++iter) {
		const SvtxTrackState* candidate = iter->second;
		float dist = sqrt(
				pow(candidate->get_x() - intersection[0], 2)
						+ pow(candidate->get_y() - intersection[1], 2)
						+ pow(candidate->get_z() - intersection[2], 2));

		if (dist < min_dist) {
			closest = candidate;
			min_dist = dist;
		}
	}

	// if we just got back the previous case, bail
	if (closest->get_pathlength() == 0.0)
		return;

	// recompute using the closer state vector
	projectToRadius(closest, track->get_charge(), B, radius, intersection);

	return;
}

void PHG4KalmanPatRec::projectToRadius(const SvtxTrackState* state, int charge,
		double B, double radius, std::vector<double>& intersection) {
	intersection.clear();
	intersection.assign(3, NAN);

	// find 2d intersections in x,y plane
	std::set<std::vector<double> > intersections;
	if (B != 0.0) {
		// magentic field present, project track as a circle leaving the state position

		// compute the center of rotation and the helix parameters
		// x(u) = r*cos(q*u+cphi) + cx
		// y(u) = r*sin(q*u+cphi) + cy
		// z(u) = b*u + posz;

		double cr = state->get_pt() * 333.6 / B;          // radius of curvature
		double cx = state->get_x()
				- (state->get_py() * cr) / charge / state->get_pt(); // center of rotation, x
		double cy = (state->get_px() * cr) / charge / state->get_pt()
				+ state->get_y(); // center of rotation, y
		double cphi = atan2(state->get_y() - cy, state->get_x() - cx); // phase of state position
		double b = state->get_pz() / state->get_pt() * cr;      // pitch along z

		if (!circle_circle_intersections(0.0, 0.0, radius, cx, cy, cr,
				&intersections)) {
			return;
		}

		if (intersections.empty())
			return;

		// downselect solutions based on track direction
		// we want the solution where the state vector would exit the cylinder
		// this can be determined by the direction that the track circulates in

		// rotate the px,py to the postion of the solution
		// then ask if the dot product of the displacement vector between the solution
		// and the cylinder center with the rotated momentum vector is positive
		std::set<std::vector<double> >::iterator remove_iter =
				intersections.end();
		double intersection_z = 0.0;
		for (std::set<std::vector<double> >::iterator iter =
				intersections.begin(); iter != intersections.end(); ++iter) {
			double x = iter->at(0);
			double y = iter->at(1);

			// find the azimuthal rotation about the center of rotation between the state vector and the solution

			// displacement between center of rotation and solution
			double dx = x - cx;
			double dy = y - cy;
			double dphi = atan2(dy, dx);

			// displacement between center of rotation and state position
			double state_dx = state->get_x() - cx;
			double state_dy = state->get_y() - cy;
			double state_dphi = atan2(state_dy, state_dx);

			// relative rotation angle
			double rotphi = (dphi - state_dphi);

			// rotated momentum at the solution
			double rotpx = cos(rotphi) * state->get_px()
					- sin(rotphi) * state->get_py();
			double rotpy = sin(rotphi) * state->get_px()
					+ cos(rotphi) * state->get_py();

			// assumes cylinder is centered at 0,0
			double dot = rotpx * x + rotpy * y;

			// our solution will have a momentum vector leaving the cylinder surface
			if (dot >= 0.0) {
				// find the solution for z
				double u = (dphi - cphi) / charge;

				// look only along the projection (not backward)
				if (u > 2.0 * M_PI) {
					u = u - 2.0 * M_PI;
				} else if (u < 0.0) {
					u = u + 2.0 * M_PI;
				}

				intersection_z = b * u + state->get_z();
			} else {
				remove_iter = iter;
			}
		}

		if (remove_iter != intersections.end()) {
			intersections.erase(remove_iter);
		}

		if (intersections.empty())
			return;

		intersection[0] = intersections.begin()->at(0);
		intersection[1] = intersections.begin()->at(1);
		intersection[2] = intersection_z;

		return;

	} else {
		// no magnetic field project track as a line

		circle_line_intersections(0.0, 0.0, radius, state->get_x(),
				state->get_y(), state->get_px(), state->get_py(),
				&intersections);

		if (intersections.empty())
			return;

		// downselect solutions based on track direction
		// we want the solution where the state vector would point outward
		// since the track doesn't bend this will be the solution where
		// the dot product of the displacement between the solution and the cylinder center
		// and the momentum vector is positive
		std::set<std::vector<double> >::iterator remove_iter =
				intersections.end();
		double intersection_z = 0.0;
		for (std::set<std::vector<double> >::iterator iter =
				intersections.begin(); iter != intersections.end(); ++iter) {
			double x = iter->at(0);
			double y = iter->at(1);

			// assumes cylinder is centered at 0,0
			double dot = state->get_px() * x + state->get_py() * y;
			if (dot >= 0.0) {

				// x(u) = px*u + x1
				// y(u) = py*u + y1
				// z(u) = pz*u + z1

				double u = NAN;
				if (state->get_px() != 0) {
					u = (intersection[0] - state->get_x()) / state->get_px();
				} else if (state->get_py() != 0) {
					u = (intersection[1] - state->get_y()) / state->get_py();
				}

				intersection_z = state->get_pz() * u + state->get_z();
			} else {
				remove_iter = iter;
			}
		}

		if (remove_iter != intersections.end()) {
			intersections.erase(remove_iter);
		}

		if (intersections.empty())
			return;

		intersection[0] = intersections.begin()->at(0);
		intersection[1] = intersections.begin()->at(1);
		intersection[2] = intersection_z;

		return;
	}

	return;
}

void PHG4KalmanPatRec::set_material(int layer, float value) {
	_user_material[layer] = value;
}

int PHG4KalmanPatRec::CreateNodes(PHCompositeNode* topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	PHNodeIterator iter_dst(dstNode);

	// Create the SVTX node
	PHCompositeNode* tb_node =
			dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
					"SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (Verbosity() > 0)
			cout << "SVTX node added" << endl;
	}

	_g4tracks = new SvtxTrackMap_v1;
	PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(_g4tracks,
			"SvtxTrackMap", "PHObject");
	tb_node->addNode(tracks_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxTrackMap node added" << endl;

	_g4vertexes = new SvtxVertexMap_v1;
	PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
			_g4vertexes, "SvtxVertexMap", "PHObject");
	tb_node->addNode(vertexes_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxVertexMap node added" << endl;

	/*
	 PHG4CylinderGeomContainer* geoms =
	 findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SVTX");
	 if (!geoms) {
	 cerr << PHWHERE << " ERROR: Can't find CYLINDERGEOM_SVTX Node." << endl;
	 return Fun4AllReturnCodes::ABORTEVENT;
	 }
	 */

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4KalmanPatRec::InitializeGeometry(PHCompositeNode *topNode) {
  
  //---------------------------------------------------------
  // Grab Run-Dependent Detector Geometry and Configure Hough
  //---------------------------------------------------------
  
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<
  PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SILICON_TRACKER");
  PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");
  
  //  if (cellgeos || laddergeos || mapsladdergeos) {
  //    unsigned int ncelllayers = 0;
  //    if (cellgeos) ncelllayers += cellgeos->get_NLayers();
  //    unsigned int nladderlayers = 0;
  //    if (laddergeos) nladderlayers += laddergeos->get_NLayers();
  //    unsigned int nmapsladderlayers = 0;
  //    if (mapsladdergeos) nmapsladderlayers += mapsladdergeos->get_NLayers();
  //    _nlayers_seeding = ncelllayers + nladderlayers + nmapsladderlayers;
  //  } else {
  //    cerr << PHWHERE
  //         << "None of  CYLINDERCELLGEOM_SVTX or CYLINDERGEOM_SILICON_TRACKER or CYLINDERGEOM_MAPS"
  //            "available, bail"
  //         << std::endl;
  //    return Fun4AllReturnCodes::ABORTRUN;
  //  }
  
  //  _nlayers_seeding = 7;
  //  int seeding_layer_array[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  //  _seeding_layer.assign(seeding_layer_array, seeding_layer_array+9 );
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
  
  // get node containing the digitized hits
  _svtxhitsmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_svtxhitsmap) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  _cells_svtx = findNode::getClass<PHG4CellContainer>(topNode,
						      "G4CELL_SVTX");
  
  _cells_intt = findNode::getClass<PHG4CellContainer>(
						      topNode, "G4CELL_SILICON_TRACKER");
  
  _cells_maps = findNode::getClass<PHG4CellContainer>(
						      topNode, "G4CELL_MAPS");
  
  if (!_cells_svtx and !_cells_intt and !_cells_maps) {
    if (Verbosity() >= 0) {
      LogError("No PHG4CellContainer found!");}
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  _geom_container_intt = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SILICON_TRACKER");
  
  _geom_container_maps = findNode::getClass<
  PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");
  
  if (!_cells_svtx && !_cells_maps && !_cells_intt) {
    cout << PHWHERE << "ERROR: Can't find any cell node!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4KalmanPatRec::InitializePHGenFit(PHCompositeNode* topNode) {

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

int PHG4KalmanPatRec::setup_seed_tracker_objects() {

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

int PHG4KalmanPatRec::setup_initial_tracker_object() {

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

int PHG4KalmanPatRec::setup_tracker_object() {

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

int PHG4KalmanPatRec::GetNodes(PHCompositeNode* topNode) {

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	// Pull the reconstructed track information off the node tree...
	_bbc_vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");

	_g4clusters = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_g4clusters) {
		cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	if(_hit_used_map_size!=0) delete[] _hit_used_map;

	_hit_used_map_size = static_cast<int>(_g4clusters->size());
	_hit_used_map = new int[_hit_used_map_size];
	for (Int_t i=0;i<_hit_used_map_size;i++){
	  _hit_used_map[i] = 0;
	}


	// Pull the reconstructed track information off the node tree...
	_g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_g4tracks) {
		cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Pull the reconstructed track information off the node tree...
	_g4vertexes = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_g4vertexes) {
		cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4KalmanPatRec::translate_input() {

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
  for (SvtxClusterMap::Iter iter = _g4clusters->begin();
       iter != _g4clusters->end(); ++iter) {
    if(_hit_used_map[iter->first]!=0){
      continue;
    }
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
      << "PHG4KalmanPatRec::process_event has the following input clusters:"
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

int PHG4KalmanPatRec::fast_vertex_from_bbc() {

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

int PHG4KalmanPatRec::fast_vertex_guessing() {

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

int PHG4KalmanPatRec::initial_vertex_finding() {

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

//FIXME this is now a simulator
int PHG4KalmanPatRec::vertexing(PHCompositeNode* topNode) {
	PHG4TruthInfoContainer* g4truth = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
	PHG4VtxPoint* first_point = g4truth->GetPrimaryVtx(
			g4truth->GetPrimaryVertexIndex());

	_vertex.clear();
	_vertex.assign(3, 0.0);

	_vertex[0] = first_point->get_x();
	_vertex[1] = first_point->get_y();
	_vertex[2] = first_point->get_z();

#ifndef __CINT__
	gsl_rng *RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned int seed = PHRandomSeed(); // fixed seed is handled in this funtcion
//  cout << Name() << " random seed: " << seed << endl;
	gsl_rng_set(RandomGenerator, seed);

//	_vertex[0] += _vertex_error[0] * gsl_ran_ugaussian(RandomGenerator);
//	_vertex[1] += _vertex_error[1] * gsl_ran_ugaussian(RandomGenerator);
//	_vertex[2] += _vertex_error[2] * gsl_ran_ugaussian(RandomGenerator);

	gsl_rng_free(RandomGenerator);
#endif

	if (Verbosity() > 1) {
		cout << __LINE__ << " PHG4KalmanPatRec::vertexing: {" << _vertex[0]
				<< ", " << _vertex[1] << ", " << _vertex[2] << "} +- {"
				<< _vertex_error[0] << ", " << _vertex_error[1] << ", "
				<< _vertex_error[2] << "}" << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4KalmanPatRec::full_track_seeding() {

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

int PHG4KalmanPatRec::export_output() {

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
			if ((track_hits.at(ihit).get_id()) >= _g4clusters->size()) {
				continue;
			}
			SvtxCluster* cluster = _g4clusters->get(
					track_hits.at(ihit).get_id());
			//mark hit asu used by iteration number n
			_hit_used_map[track_hits.at(ihit).get_id()] = _n_iteration;
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

		_g4tracks->insert(&track);
		vertex.insert_track(track.get_id());

		if (Verbosity() > 5) {
			cout << "track " << itrack << " quality = " << track.get_quality()
					<< endl;
			cout << "px = " << track.get_px() << " py = " << track.get_py()
					<< " pz = " << track.get_pz() << endl;
		}
	}  // track loop

	SvtxVertex *vtxptr = _g4vertexes->insert(&vertex);
	if (Verbosity() > 5)
		vtxptr->identify();

	if (Verbosity() > 0) {
		cout << "PHG4KalmanPatRec::process_event -- leaving process_event"
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

int PHG4KalmanPatRec::add_tracks() {

	if (_tracks.empty())
		return Fun4AllReturnCodes::EVENT_OK;

	vector<SimpleHit3D> track_hits;

	//Mark used clusters
	for (unsigned int itrack = 0; itrack < _tracks.size(); itrack++) {
	  _all_tracks.push_back(_tracks[itrack]);
	  _all_track_errors.push_back(_track_errors[itrack]);
	  _all_track_covars.push_back(_track_covars[itrack]);
	  
	  track_hits = _tracks.at(itrack).hits;
	  
	  for (unsigned int ihit = 0; ihit < track_hits.size(); ihit++) {
	    if ((track_hits.at(ihit).get_id()) >= _g4clusters->size()) {
	      continue;
	    }
	    //mark hit as used by iteration number n
	    _hit_used_map[track_hits.at(ihit).get_id()] = _n_iteration;
	  }

	}  // track loop

	// we are done with these now...
	_clusters.clear();
	_tracks.clear();
	_track_errors.clear();
	_track_covars.clear();

	return Fun4AllReturnCodes::EVENT_OK;
}

float PHG4KalmanPatRec::kappaToPt(float kappa) {
	return _pt_rescale * _magField / 333.6 / kappa;
}

float PHG4KalmanPatRec::ptToKappa(float pt) {
	return _pt_rescale * _magField / 333.6 / pt;
}

void PHG4KalmanPatRec::convertHelixCovarianceToEuclideanCovariance(float B,
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

void PHG4KalmanPatRec::shift_coordinate_system(double dx, double dy,
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

bool PHG4KalmanPatRec::circle_line_intersections(double x0, double y0,
		double r0, double x1, double y1, double vx1, double vy1,
		std::set<std::vector<double> >* points) {
	// P0: center of rotation
	// P1: point on line
	// P2: second point on line
	// P3: intersections

	// dr: distance between P1 & P2
	// delta: discriminant on number of solutions

	points->clear();

	double x2 = x1 + vx1;
	double y2 = y1 + vy1;

	double dr = sqrt(pow(vx1, 2) + pow(vy1, 2));
	double det = x1 * y2 - x2 * y1;

	double delta = pow(r0, 2) * pow(dr, 2) - pow(det, 2);
	if (delta < 0)
		return false;

	double sgn_vy1 = 1.0;
	if (vy1 < 0.0)
		sgn_vy1 = -1.0;

	double x3 = (det * vy1
			+ sgn_vy1 * vx1 * sqrt(pow(r0, 2) * pow(dr, 2) - pow(det, 2)))
			/ pow(dr, 2);
	double y3 = (-1.0 * det * vx1
			+ fabs(vy1) * sqrt(pow(r0, 2) * pow(dr, 2) - pow(det, 2)))
			/ pow(dr, 2);

	std::vector<double> p3;
	p3.push_back(x3);
	p3.push_back(y3);
	points->insert(p3);

	x3 = (det * vy1
			- sgn_vy1 * vx1 * sqrt(pow(r0, 2) * pow(dr, 2) - pow(det, 2)))
			/ pow(dr, 2);
	y3 = (-1.0 * det * vx1
			- fabs(vy1) * sqrt(pow(r0, 2) * pow(dr, 2) - pow(det, 2)))
			/ pow(dr, 2);

	p3[0] = x3;
	p3[1] = y3;
	points->insert(p3);

	return true;
}

int PHG4KalmanPatRec::CleanupSeedsByHitPattern() {

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

int PHG4KalmanPatRec::CleanupTracksByHitPattern() {

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


int PHG4KalmanPatRec::check_track_exists(MapPHGenFitTrack::iterator iter){
	
  
  //Loop over hitIDs on current track and check if they have been used
  unsigned int n_clu = iter->second->get_cluster_IDs().size();

  unsigned int  n_clu_used = 0;
  const std::vector<unsigned int>& clusterIDs = iter->second->get_cluster_IDs();
  for(unsigned int iCluId = 0; iCluId < clusterIDs.size(); ++iCluId){
    unsigned int cluster_ID = clusterIDs[iCluId];
    if(_hit_used_map[cluster_ID]>0)n_clu_used++;
  }
  int code = 0;
  if(((float)n_clu_used/n_clu)>0.3){
    if(Verbosity()>=1)
      cout << "Found duplicate track. n_clu: " << n_clu << " c_clu_used: " << n_clu_used << " n_iter: " << _n_iteration<< endl;
    /*
    for(unsigned int iCluId = 0; iCluId < clusterIDs.size(); ++iCluId){
      unsigned int cluster_ID = clusterIDs[iCluId];
      cout << "#Clu_g = " << iCluId 
	   << " layer: " << _g4clusters->get(cluster_ID)->get_layer() 
	   << " r: " << TMath::Sqrt(_g4clusters->get(cluster_ID)->get_x()*_g4clusters->get(cluster_ID)->get_x() +_g4clusters->get(cluster_ID)->get_y()*_g4clusters->get(cluster_ID)->get_y() )
	   << endl;
    }
    */
    return code;
  }
  code = 1;
  return code;
}

int PHG4KalmanPatRec::CleanupSeeds() {

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

int PHG4KalmanPatRec::FullTrackFitting(PHCompositeNode* topNode) {

#ifdef _DEBUG_
	std::cout << "=========================" << std::endl;
	std::cout << "PHG4KalmanPatRec::FullTrackFitting: Start: Event: "<< _event << std::endl;
	std::cout << "Total Raw Tracks: " << _tracks.size() << std::endl;
	std::cout << "=========================" << std::endl;
#endif

	/*!
	 *   sort clusters
	 */
	BuildLayerZPhiHitMap();

	vector<genfit::Track*> evt_disp_copy;

	_PHGenFitTracks.clear();

	//	for(unsigned int itrack = 0; itrack < ( (_tracks.size() < 100) ? _tracks.size() : 100); ++itrack) {
	for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
#ifdef _DEBUG_
		std::cout
		<< __LINE__
		<< ": Processing itrack: " << itrack
		<< ": Total tracks: " << _g4tracks->size()
		<<endl;
#endif

		/*!
		 * Translate SimpleTrack3D To PHGenFitTracks
		 */
		if(Verbosity() > 1) _t_translate_to_PHGenFitTrack->restart();
		SimpleTrack3DToPHGenFitTracks(topNode, itrack);
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
			<< ": propergating: " << i <<"/" << itrack
			<< endl;
#endif

			std::vector<unsigned int> clusterIDs = iter->second->get_cluster_IDs();

			unsigned int init_layer = UINT_MAX;

			if(!is_splitting_track) {
				if(_init_direction == 1) {
					init_layer = _g4clusters->get(clusterIDs.front())->get_layer();
					TrackPropPatRec(topNode, iter, init_layer, _nlayers_all, true);
					TrackPropPatRec(topNode, iter, init_layer, 0, false);
				} else {
					init_layer = _g4clusters->get(clusterIDs.back())->get_layer();
					TrackPropPatRec(topNode, iter, init_layer, 0, true);
					TrackPropPatRec(topNode, iter, init_layer, _nlayers_all, false);
				}
				is_splitting_track = true;
			} else {
				if(_init_direction == 1) {
					init_layer = _g4clusters->get(clusterIDs.front())->get_layer();
					TrackPropPatRec(topNode, iter, init_layer, _nlayers_all, false);
				} else {
					init_layer = _g4clusters->get(clusterIDs.back())->get_layer();
					TrackPropPatRec(topNode, iter, init_layer, 0, false);
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
			OutputPHGenFitTrack(topNode, iter);
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
	std::cout << "PHG4KalmanPatRec::FullTrackFitting: End: Event: "<< _event << std::endl;
	std::cout << "Total Final Tracks: " << _g4tracks->size() << std::endl;
	std::cout << "=========================" << std::endl;
#endif

	if (_do_evt_display) {
		_fitter->getEventDisplay()->addEvent(evt_disp_copy);
	} else {
		evt_disp_copy.clear();
	}

	_tracks.clear();

	return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4KalmanPatRec::ExportOutput() { return 0;}

int PHG4KalmanPatRec::OutputPHGenFitTrack(PHCompositeNode* topNode, MapPHGenFitTrack::iterator iter) {
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
//#ifdef _DEBUG_
//	std::cout << "=========================" << std::endl;
//	std::cout << "PHG4KalmanPatRec::FullTrackFitting: Event: "<< _event << std::endl;
//	std::cout << "Total Raw Tracks: " << _trackID_PHGenFitTrack.size() << std::endl;
//	std::cout << "=========================" << std::endl;
//#endif

//	for (std::map<int, std::shared_ptr<PHGenFit::Track>>::iterator iter =
//			_trackID_PHGenFitTrack.begin();
//			iter != _trackID_PHGenFitTrack.end(); iter++) {

#ifdef _DEBUG_
		std::cout << "=========================" << std::endl;
		//std::cout << __LINE__ << ": iPHGenFitTrack: " << iter->first << std::endl;
		std::cout << __LINE__ << ": _g4tracks->size(): " << _g4tracks->size() << std::endl;
		std::cout << "Contains: " << iter->second->get_cluster_IDs().size() << " clusters." <<std::endl;
		std::cout << "=========================" << std::endl;
#endif

		SvtxTrack_v1 track;
		//track.set_id(iter->first);
		track.set_id(_g4tracks->size());

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
		  SvtxCluster* cluster = clustermap->get(cluster_id);
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
		//		if(is_good_track||_n_iteration==4)
		//if(is_good_track||_n_iteration>=0)
		if(_n_iteration>=0)
		  {
		    for(unsigned int cluster_ID : iter->second->get_cluster_IDs()){
		      _hit_used_map[cluster_ID] = _n_iteration;
		    }
		    
		    _g4tracks->insert(&track);
		  }
		if (Verbosity() > 5) {
			cout << "track " << _g4tracks->size() << " quality = " << track.get_quality()
					<< endl;
			cout << "px = " << track.get_px() << " py = " << track.get_py()
					<< " pz = " << track.get_pz() << endl;
		}

		if(Verbosity() >= 1) _t_output_io->stop();
//	}



	return 0;
}


int PHG4KalmanPatRec::SimpleTrack3DToPHGenFitTracks(PHCompositeNode* topNode, unsigned int itrack) {

	// clean up working array for each event
	_PHGenFitTracks.clear();

	double time1 = 0;
	double time2 = 0;

	if(Verbosity() > 1){
	  time1 = _t_translate1->get_accumulated_time();
	  time2 = _t_translate1->get_accumulated_time();
	  _t_translate1->restart();
	}

	vector<SimpleHit3D> track_hits = _tracks.at(itrack).hits;

	float kappa = _tracks.at(itrack).kappa;
	float d = _tracks.at(itrack).d;
	float phi = _tracks.at(itrack).phi;
	float dzdl = _tracks.at(itrack).dzdl;
	float z0 = _tracks.at(itrack).z0;
	float nhit = (float) _tracks.at(itrack).hits.size();
	float ml = 0;
	float rec = 0;
	float dt = 0;

	if(_tracks.at(itrack).hits.size() == 0) {
	  if(Verbosity() > 1) _t_translate1->stop();
	  return 1;
	}

	if(!(kappa==kappa && d==d && phi==phi && dzdl==dzdl && z0==z0)) {
	  if(Verbosity() > 1) _t_translate1->stop();
	  return 1;
	}

	if(kappa == 0) {
	  if(Verbosity() > 1) _t_translate1->stop();
	  return 1;
	}

	float pT = kappaToPt(kappa);


	//FIXME
	if(pT < _cut_min_pT) return 1;

	float x_center = cos(phi) * (d + 1 / kappa); // x coordinate of circle center
	float y_center = sin(phi) * (d + 1 / kappa);  // y    "      "     " "


	// find helicity from cross product sign
	short int helicity;
	{
		unsigned int hitID0 = track_hits.front().get_id();
		unsigned int hitID1 = track_hits.back().get_id();
		SvtxCluster* cluster0 = _g4clusters->get(hitID0);
		SvtxCluster* cluster1 = _g4clusters->get(hitID1);

		if ((cluster0->get_x() - x_center) * (cluster1->get_y() - y_center)
				- (cluster0->get_y() - y_center)
						* (cluster1->get_x() - x_center) > 0) {
			helicity = 1;
		} else {
			helicity = -1;
		}
	}

	float pZ = 0;
	if (dzdl != 1) {
		pZ = pT * dzdl / sqrt(1.0 - dzdl * dzdl);
	}

	TVector3 seed_mom(
			pT * cos(phi - helicity * M_PI / 2),
			pT * sin(phi - helicity * M_PI / 2),
			pZ
	);

	TVector3 seed_pos(
			_vertex[0] + d * cos(phi),
			_vertex[1] + d * sin(phi),
			_vertex[2] + z0
			);

	Eigen::Matrix<float, 6, 6> euclidean_cov =
			Eigen::Matrix<float, 6, 6>::Zero(6, 6);
	convertHelixCovarianceToEuclideanCovariance(_magField, phi, d, kappa,
			z0, dzdl, _track_covars[itrack], euclidean_cov);

	if(Verbosity() > 1) _t_translate2->restart();

	//TODO optimize
	const float blowup_factor = 1.;

	TMatrixDSym seed_cov(6);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			seed_cov[i][j] = blowup_factor*euclidean_cov(i, j);
		}
	}

	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
	std::shared_ptr<PHGenFit::Track> track(
			new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov));

	std::multimap<float, unsigned int> m_r_clusterID;
	for (SimpleHit3D hit : track_hits) {

		unsigned int cluster_ID = hit.get_id();
		SvtxCluster* cluster = _g4clusters->get(cluster_ID);

		float r = sqrt(
				cluster->get_x() * cluster->get_x() +
				cluster->get_y() * cluster->get_y());

		m_r_clusterID.insert(std::pair<float, unsigned int>(r, hit.get_id()));
	}

	std::vector<PHGenFit::Measurement*> measurements;
	{
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

		SvtxCluster* cluster = _g4clusters->get(cluster_ID);
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
		  _analyzing_ntuple->Fill(pT,kappa,d,phi,dzdl,z0,nhit,ml/nhit,rec,dt);
		return -1;
	}

	int nhits = track->get_cluster_IDs().size();
	float chi2 = track->get_chi2();
	float ndf  = track->get_ndf();

	if(nhits > 0 and chi2 > 0 and ndf > 0) {
		_PHGenFitTracks.push_back(
				MapPHGenFitTrack::value_type(
						PHG4KalmanPatRec::TrackQuality(nhits, chi2, ndf, nhits, 0, 0), track)
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
	  _analyzing_ntuple->Fill(pT,kappa,d,phi,dzdl,z0,nhit,rec,dt);

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4KalmanPatRec::TrackPropPatRec(
		PHCompositeNode* topNode,
		//const int iPHGenFitTrack, std::shared_ptr<PHGenFit::Track> &track,
		MapPHGenFitTrack::iterator &track_iter,
		unsigned int init_layer, unsigned int end_layer,
		const bool use_fitted_state_once) {

	std::shared_ptr<PHGenFit::Track> &track = track_iter->second;

        //some debug info
	/*
	  float init_pt = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom().Pt();
	  float init_phi = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom().Phi();
	  float init_eta = 0.0;//track->getGenFitTrack()->getCardinalRep()->getMom().Eta();
	*/
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
			if (_g4clusters->get(clusterIDs[i])->get_layer() == init_layer) {
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
				SvtxCluster* extrapolate_base_cluster = _g4clusters->get(extrapolate_base_cluster_id);
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
			SvtxCluster* cluster = _g4clusters->get(cluster_ID);
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
		
		PHG4KalmanPatRec::TrackQuality tq(track_iter->first);

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
								PHG4KalmanPatRec::TrackQuality(
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

PHGenFit::Measurement* PHG4KalmanPatRec::SvtxClusterToPHGenFitMeasurement(
		const SvtxCluster* cluster) {

	if(!cluster) return nullptr;

	TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
	TVector3 n(cluster->get_x(), cluster->get_y(), 0);

	unsigned int begin_hit_id = *(cluster->begin_hits());
	//LogDebug(begin_hit_id);
	SvtxHit* svtxhit = _svtxhitsmap->find(begin_hit_id)->second;
	//LogDebug(svtxhit->get_cellid());

	PHG4Cell* cell_svtx = nullptr;
	PHG4Cell* cell_intt = nullptr;
	PHG4Cell* cell_maps = nullptr;

	if(_cells_svtx) cell_svtx = _cells_svtx->findCell(svtxhit->get_cellid());
	if(_cells_intt) cell_intt = _cells_intt->findCell(svtxhit->get_cellid());
	if(_cells_maps) cell_maps = _cells_maps->findCell(svtxhit->get_cellid());
	if(!(cell_svtx or cell_intt or cell_maps)){
		if(Verbosity()>=0)
			LogError("!(cell_svtx or cell_intt or cell_maps)");
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
		PHG4CylinderGeom_MAPS *geom =
				(PHG4CylinderGeom_MAPS*) _geom_container_maps->GetLayerGeom(
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
		PHG4CylinderGeomSiLadders* geom =
		  dynamic_cast<PHG4CylinderGeomSiLadders*> (_geom_container_intt->GetLayerGeom(
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

int PHG4KalmanPatRec::BuildLayerZPhiHitMap() {

	_layer_thetaID_phiID_cluserID.clear();

	//for(SimpleHit3D cluster : _clusters){
	for (SvtxClusterMap::Iter iter = _g4clusters->begin();
			iter != _g4clusters->end(); ++iter) {
		SvtxCluster* cluster = iter->second;

		unsigned int layer = cluster->get_layer();

		float x = cluster->get_x()-_vertex[0];
		float y = cluster->get_y()-_vertex[1];
		float z = cluster->get_z()-_vertex[2];

		float phi = atan2(y,x);
		float r = sqrt(x*x + y*y);
		float theta = atan2(r,z);

#ifdef _DEBUG_
		//float rphi = r*phi;
		std::cout
		<< __LINE__
		<<": ID: " << cluster->get_id()
		<<": layer: "<<cluster->get_layer()
		<<", r: "<<r
		<<", z: "<<z
		<<", phi: "<<phi
		<<", theta: "<<theta
		<<endl;
#endif

		unsigned int idx = encode_cluster_index(layer, theta, phi);

#ifdef _DEBUG_
			cout
			<<__LINE__<<": "
			<<"{ "
			<<layer <<", "
			<<z <<", "
			<<phi << "} =>"
			<<idx << ": size: "
			<<_layer_thetaID_phiID_cluserID.count(idx)
			<<endl;
#endif

		_layer_thetaID_phiID_cluserID.insert(std::make_pair(idx, cluster->get_id()));
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<unsigned int> PHG4KalmanPatRec::SearchHitsNearBy(const unsigned int layer,
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
				SvtxCluster* cluster = _g4clusters->get(iter->second);
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

//std::shared_ptr<SvtxTrack> PHG4KalmanPatRec::MakeSvtxTrack(
//		const int genfit_track_ID, const SvtxVertex* vertex) {
//
//	std::shared_ptr<SvtxTrack> svtxtrack(new SvtxTrack_v1());
//
//	return svtxtrack;
//}

unsigned int PHG4KalmanPatRec::encode_cluster_index(const unsigned int layer,
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
	std::cout<<__LINE__<<": "
			<<": layer: "<<layer
			<<", irphi: "<<irphi
			<<", itheta: "<<itheta
			<<endl;
#endif

	return encode_cluster_index(layer, itheta, irphi);
}

//unsigned int PHG4KalmanPatRec::encode_cluster_index(const unsigned int layer,
//		const unsigned int iz, const unsigned int irphi) {
//
//	std::bitset<7> layer_bits(layer);
//	std::bitset<11> z_bits(iz);
//	std::bitset<14> rphi_bits(irphi);
//	std::bitset<32> idx_bits(0);
//
//	for(unsigned int i=0;i<idx_bits.size();++i){
//		if(i < 14)
//			idx_bits[i] = rphi_bits[i];
//		else if(i < 25)
//			idx_bits[i] = z_bits[i-14];
//		else
//			idx_bits[i] = layer_bits[i-25];
//	}
//
//	return (unsigned int) idx_bits.to_ulong();
//}

unsigned int PHG4KalmanPatRec::encode_cluster_index(const unsigned int layer,
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

bool PHG4KalmanPatRec::circle_circle_intersections(double x0, double y0,
		double r0, double x1, double y1, double r1,
		std::set<std::vector<double> >* points) {
	// P0: center of rotation on first circle
	// P1: center of rotation of second circle
	// P2: point between P0 & P1 on radical line
	// P3: intersection points

	// d: distance between P0 and P1
	// a: distance between P0 and P2
	// h: distance from P2 and P3s

	points->clear();

	// distance between two circle centers
	double d = sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));

	// handle error conditions
	if (fabs(r0 + r1) < d)
		return false; // no solution
	if (fabs(r0 - r1) > d)
		return false; // no solution
	if (d == 0 && r0 == r1)
		return false; // infinite solutions

	// compute distances to intersection points
	double a = (pow(r0, 2) - pow(r1, 2) + pow(d, 2)) / (2 * d);
	double h = sqrt(pow(r0, 2) - pow(a, 2));

	// compute P2
	double x2 = x0 + a * (x1 - x0) / d;
	double y2 = y0 + a * (y1 - y0) / d;

	// compute intersection, p3
	double x3 = x2 + h * (y1 - y0) / d;
	double y3 = y2 - h * (x1 - x0) / d;

	std::vector<double> p3;
	p3.push_back(x3);
	p3.push_back(y3);
	points->insert(p3);

	// second intersection (if different than first)
	x3 = x2 - h * (y1 - y0) / d;
	y3 = y2 + h * (x1 - x0) / d;

	p3[0] = x3;
	p3[1] = y3;
	points->insert(p3);

	return points;
}


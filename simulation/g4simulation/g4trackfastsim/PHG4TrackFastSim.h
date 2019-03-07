/*!
 *  \file		PHG4TrackFastSim.h
 *  \brief		Kalman Filter based on smeared truth PHG4Hit
 *  \details	Kalman Filter based on smeared truth PHG4Hit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4TrackFastSim_H__
#define __PHG4TrackFastSim_H__

#include <iostream>
#include <string>
#include <vector>
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>
#include <phgenfit/Measurement.h>

#include <TMatrixDSym.h>
#include <TVector3.h>

class PHG4Particle;
namespace PHGenFit {
class PlanarMeasurement;
} /* namespace PHGenFit */

namespace PHGenFit {
class Track;
} /* namespace PHGenFit */

namespace genfit {
class GFRaveVertexFactory;
} /* namespace genfit */

class SvtxTrack;
namespace PHGenFit {
class Fitter;
} /* namespace PHGenFit */

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;

class PHG4TrackFastSim: public SubsysReco {
public:

	enum DETECTOR_TYPE {Vertical_Plane, Cylinder};

	//! Default constructor
	PHG4TrackFastSim(const std::string &name = "PHG4TrackFastSim");

	//! dtor
	~PHG4TrackFastSim();

	//!Initialization, called for initialization
	int Init(PHCompositeNode *);

	//!Initialization Run, called for initialization of a run
	int InitRun(PHCompositeNode *);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *);

	//!End, write and close files
	int End(PHCompositeNode *);

	bool is_do_evt_display() const {
		return _do_evt_display;
	}

	void set_do_evt_display(bool doEvtDisplay) {
		_do_evt_display = doEvtDisplay;
	}

	double get_phi_resolution() const {
		return _phi_resolution;
	}

	void set_phi_resolution(double fgemPhiResolution) {
		_phi_resolution = fgemPhiResolution;
	}

	double get_r_resolution() const {
		return _r_resolution;
	}

	void set_r_resolution(double fgemRResolution) {
		_r_resolution = fgemRResolution;
	}

	const std::string& get_fit_alg_name() const {
		return _fit_alg_name;
	}

	void set_fit_alg_name(const std::string& fitAlgName) {
		_fit_alg_name = fitAlgName;
	}

	double get_pat_rec_hit_finding_eff() const {
		return _pat_rec_hit_finding_eff;
	}

	void set_pat_rec_hit_finding_eff(double patRecHitFindingEff) {
		if(!(patRecHitFindingEff>=0&&patRecHitFindingEff<=1)) {
			std::cout<<"ERROR: _pat_rec_hit_finding_eff out of range! \n";
		}
		_pat_rec_hit_finding_eff = patRecHitFindingEff;
	}

	double get_pat_rec_noise_prob() const {
		return _pat_rec_noise_prob;
	}

	void set_pat_rec_noise_prob(double patRecNoiseProb) {
		if(!(patRecNoiseProb <= 1. && patRecNoiseProb >= 0)) {
			std::cout<<"ERROR: _pat_rec_noise_prob out of range! \n";
			return;
		}
		_pat_rec_noise_prob = patRecNoiseProb;
	}

	const std::vector<std::string>& get_phg4hits_names() const {
		return _phg4hits_names;
	}

	void set_phg4hits_names(const std::vector<std::string>& phg4hitsNames) {
	  _phg4hits_names.clear();
	  _phg4_detector_type.clear();
	  _phg4_detector_radres.clear();
	  _phg4_detector_phires.clear();
	  _phg4_detector_lonres.clear();
	  _phg4_detector_hitfindeff.clear();
	  _phg4_detector_noise.clear();
	  //==
	  _phg4hits_names = phg4hitsNames;
	  _N_DETECTOR_LAYER = _phg4hits_names.size();
	  for(int i=0; i!=_N_DETECTOR_LAYER; ++i) { //defaulting
	    _phg4_detector_type.push_back(Cylinder);
	    _phg4_detector_radres.push_back(1e-4);
	    _phg4_detector_phires.push_back(1e-4);
	    _phg4_detector_lonres.push_back(1e-4);
	    _phg4_detector_hitfindeff.push_back(1.0);
	    _phg4_detector_noise.push_back(0.0);
	  }
	}

	void set_phg4hits_names(const std::string* phg4hitsNames, const int nlayer) {
	  _phg4hits_names.clear();
	  _phg4_detector_type.clear();
	  _phg4_detector_radres.clear();
	  _phg4_detector_phires.clear();
	  _phg4_detector_lonres.clear();
	  _phg4_detector_hitfindeff.clear();
	  _phg4_detector_noise.clear();
	  for(int i=0;i<nlayer;++i) {
	    _phg4hits_names.push_back(phg4hitsNames[i]);
	    _phg4_detector_type.push_back(Cylinder);
	    _phg4_detector_radres.push_back(1e-4);
	    _phg4_detector_phires.push_back(1e-4);
	    _phg4_detector_lonres.push_back(1e-4);
	    _phg4_detector_hitfindeff.push_back(1.0);
	    _phg4_detector_noise.push_back(0.0);
	  }
	  _N_DETECTOR_LAYER = _phg4hits_names.size();
	}

	void set_phg4hits_names(const std::string* phg4hitsNames,
				const DETECTOR_TYPE* phg4dettype, 
				const float *radres,
				const float *phires,
				const float *lonres,
				const float *eff,
				const float *noise,
				const int nlayer) {
	  _phg4hits_names.clear();
	  _phg4_detector_type.clear();
	  _phg4_detector_radres.clear();
	  _phg4_detector_phires.clear();
	  _phg4_detector_lonres.clear();
	  _phg4_detector_hitfindeff.clear();
	  _phg4_detector_noise.clear();
	  for(int i=0;i<nlayer;++i) {
	    _phg4hits_names.push_back(phg4hitsNames[i]);
	    _phg4_detector_type.push_back(phg4dettype[i]);
            _phg4_detector_radres.push_back(radres[i]);
            _phg4_detector_phires.push_back(phires[i]);
            _phg4_detector_lonres.push_back(lonres[i]);
	    _phg4_detector_hitfindeff.push_back(eff[i]);
	    _phg4_detector_noise.push_back(noise[i]);
	  }
	  _N_DETECTOR_LAYER = _phg4hits_names.size();
	}

	void set_state_names(const std::string* stateNames, const int nlayer) {
		_state_names.clear();
		for(int i=0;i<nlayer;++i) {
			_state_names.push_back(stateNames[i]);
		}
		_N_STATES = _state_names.size();
	}

	double get_z_resolution() const {
		return _z_resolution;
	}

	void set_z_resolution(double resolution) {
		_z_resolution = resolution;
	}

	DETECTOR_TYPE get_detector_type() const {
		return _detector_type;
	}

	void set_detector_type(DETECTOR_TYPE detectorType) {
		_detector_type = detectorType;
	}

	const std::string& get_trackmap_out_name() const {
		return _trackmap_out_name;
	}

	void set_trackmap_out_name(const std::string& trackmapOutName) {
		_trackmap_out_name = trackmapOutName;
	}

	const std::string& get_sub_top_node_name() const {
		return _sub_top_node_name;
	}

	void set_sub_top_node_name(const std::string& subTopNodeName) {
		_sub_top_node_name = subTopNodeName;
	}

	bool is_use_vertex_in_fitting() const {
		return _use_vertex_in_fitting;
	}

	void set_use_vertex_in_fitting(bool useVertexInFitting) {
		_use_vertex_in_fitting = useVertexInFitting;
	}

	double get_vertex_xy_resolution() const {
		return _vertex_xy_resolution;
	}

	void set_vertex_xy_resolution(double vertexXyResolution) {
		_vertex_xy_resolution = vertexXyResolution;
	}

	double get_vertex_z_resolution() const {
		return _vertex_z_resolution;
	}

	void set_vertex_z_resolution(double vertexZResolution) {
		_vertex_z_resolution = vertexZResolution;
	}

	int get_primary_assumption_pid() const {
		return _primary_assumption_pid;
	}

	void set_primary_assumption_pid(int primaryAssumptionPid) {
		_primary_assumption_pid = primaryAssumptionPid;
	}

	void set_primary_tracking(int pTrk) {
	        _primary_tracking = pTrk; 
	}

private:

	/*!
	 * Create needed nodes.
	 */
	int CreateNodes(PHCompositeNode *);

	/*!
	 * Get all the all the required nodes off the node tree.
	 */
	int GetNodes(PHCompositeNode *);

	/*!
	 *
	 */
	int PseudoPatternRecognition(const PHG4Particle* particle,
			std::vector<PHGenFit::Measurement*> & meas_out, TVector3& seed_pos,
			TVector3& seed_mom, TMatrixDSym& seed_cov, const bool do_smearing = true);

	PHGenFit::PlanarMeasurement* PHG4HitToMeasurementVerticalPlane(const PHG4Hit* g4hit, const double phi_resolution, const double r_resolution);

	PHGenFit::PlanarMeasurement* PHG4HitToMeasurementCylinder(const PHG4Hit* g4hit, const double phi_resolution, const double z_resolution);

	PHGenFit::Measurement* VertexMeasurement(const TVector3 &vtx, double dxy, double dz);

	/*!
	 * Make SvtxTrack from PHGenFit::Track
	 */
	SvtxTrack* MakeSvtxTrack(const PHGenFit::Track* phgf_track_in, 
				 const unsigned int truth_track_id = UINT_MAX,
				 const unsigned int nmeas = 0, const TVector3 &vtx = TVector3(0.0,0.0,0.0));

	//! Event counter
	int _event;

	DETECTOR_TYPE _detector_type; // deprecated

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;

	std::vector<PHG4HitContainer*> _phg4hits;
	std::vector<std::string> _phg4hits_names;
	std::vector<DETECTOR_TYPE> _phg4_detector_type;
	std::vector<float> _phg4_detector_radres;
	std::vector<float> _phg4_detector_phires;
	std::vector<float> _phg4_detector_lonres;
	std::vector<float> _phg4_detector_hitfindeff;
	std::vector<float> _phg4_detector_noise;

	//! Output Node pointers

	std::string _sub_top_node_name;
//	std::string _clustermap_out_name;
	std::string _trackmap_out_name;

//	SvtxClusterMap* _clustermap_out;

	SvtxTrackMap* _trackmap_out;



	/*!
	 *	GenFit fitter interface
	 */
	PHGenFit::Fitter* _fitter;


	/*!
	 * Available choices:
	 * KalmanFitter
	 * KalmanFitterRefTrack
	 * DafSimple
	 * DafRef
	 */
	std::string _fit_alg_name;

	//!
	int _primary_assumption_pid;

	//!
	bool _do_evt_display;

	/*!
	 * For PseudoPatternRecognition function.
	 */

	bool _use_vertex_in_fitting;

	double _vertex_xy_resolution;
	double _vertex_z_resolution;

	double _phi_resolution; //deprecated

	double _r_resolution; //deprecated

	double _z_resolution; //deprecated

	//!
	double _pat_rec_hit_finding_eff; //deprecated

	//!
	double _pat_rec_noise_prob; //deprecated

	//!
	int _N_DETECTOR_LAYER;

	//! 
	int _primary_tracking; 

	//!
	int _N_STATES;
	std::vector<std::string> _state_names;
	std::vector<double> _state_location;


};

#endif /*__PHG4TrackFastSim_H__*/






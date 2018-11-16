/*!
 *  \file		PHG4TrackKalmanFitter.h
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4TrackKalmanFitter_H__
#define __PHG4TrackKalmanFitter_H__

#include <fun4all/SubsysReco.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/Track.h>
#include <string>
#include <vector>

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

//! \brief Helper class for using RAVE vertex finder.
class PHRaveVertexFactory;

//! \brief		Refit SvtxTracks with PHGenFit.
class PHG4TrackKalmanFitter: public SubsysReco {
public:

	/*!
	 * OverwriteOriginalNode: default mode, overwrite original node
	 * MakeNewNode: Output extra new refit nodes
	 * DebugMode: overwrite original node also make extra new refit nodes
	 */
	enum OutPutMode {MakeNewNode, OverwriteOriginalNode, DebugMode};

	enum DetectorType {MIE, MAPS_TPC, MAPS_IT_TPC, LADDER_MAPS_TPC, LADDER_MAPS_IT_TPC, LADDER_MAPS_LADDER_IT_TPC, MAPS_LADDER_IT_TPC};

	//! Default constructor
	PHG4TrackKalmanFitter(const std::string &name = "PHG4TrackKalmanFitter");

	//! dtor
	~PHG4TrackKalmanFitter();

	//!Initialization, called for initialization
	int Init(PHCompositeNode *);

	//!Initialization Run, called for initialization of a run
	int InitRun(PHCompositeNode *);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *);

	//!End, write and close files
	int End(PHCompositeNode *);

	//Flags of different kinds of outputs
	enum Flag {
		//all disabled
		NONE = 0,
	};

	//Set the flag
	//Flags should be set like set_flag(PHG4TrackKalmanFitter::TRUTH, true) from macro
	void set_flag(const Flag& flag, const bool& value) {
		if (value)
			_flags |= flag;
		else
			_flags &= (~flag);
	}

	//! For evalution
	//! Change eval output filename
	void set_eval_filename(const char* file) {
		if (file)
			_eval_outname = file;
	}
	std::string get_eval_filename() const {
			return _eval_outname;
	}

	void fill_eval_tree(PHCompositeNode*);
	void init_eval_tree();
	void reset_eval_variables();

	bool is_do_eval() const {
		return _do_eval;
	}

	void set_do_eval(bool doEval) {
		_do_eval = doEval;
	}

	bool is_do_evt_display() const {
		return _do_evt_display;
	}

	void set_do_evt_display(bool doEvtDisplay) {
		_do_evt_display = doEvtDisplay;
	}

	const std::string& get_vertexing_method() const {
		return _vertexing_method;
	}

	void set_vertexing_method(const std::string& vertexingMethod) {
		_vertexing_method = vertexingMethod;
	}

	bool is_fit_primary_tracks() const {
		return _fit_primary_tracks;
	}

	void set_fit_primary_tracks(bool fitPrimaryTracks) {
		_fit_primary_tracks = fitPrimaryTracks;
	}

	OutPutMode get_output_mode() const {
		return _output_mode;
	}

	/*!
	 * set output mode, default is OverwriteOriginalNode
	 */
	void set_output_mode(OutPutMode outputMode) {
		_output_mode = outputMode;
	}

	const std::string& get_track_fitting_alg_name() const {
		return _track_fitting_alg_name;
	}

	void set_track_fitting_alg_name(const std::string& trackFittingAlgName) {
		_track_fitting_alg_name = trackFittingAlgName;
	}

	int get_primary_pid_guess() const {
		return _primary_pid_guess;
	}

	void set_primary_pid_guess(int primaryPidGuess) {
		_primary_pid_guess = primaryPidGuess;
	}

	double get_fit_min_pT() const {
		return _fit_min_pT;
	}

	void set_fit_min_pT(double cutMinPT) {
		_fit_min_pT = cutMinPT;
	}

	bool is_over_write_svtxtrackmap() const {
		return _over_write_svtxtrackmap;
	}

	void set_over_write_svtxtrackmap(bool overWriteSvtxtrackmap) {
		_over_write_svtxtrackmap = overWriteSvtxtrackmap;
	}

	bool is_over_write_svtxvertexmap() const {
		return _over_write_svtxvertexmap;
	}

	void set_over_write_svtxvertexmap(bool overWriteSvtxvertexmap) {
		_over_write_svtxvertexmap = overWriteSvtxvertexmap;
	}

	bool is_use_truth_vertex() const {
		return _use_truth_vertex;
	}

	void set_use_truth_vertex(bool useTruthVertex) {
		_use_truth_vertex = useTruthVertex;
	}

	double get_vertex_min_ndf() const {
		return _vertex_min_ndf;
	}

	void set_vertex_min_ndf(double vertexMinPT) {
		_vertex_min_ndf = vertexMinPT;
	}

private:

	//! Event counter
	int _event;

	//! Get all the nodes
	int GetNodes(PHCompositeNode *);

	//!Create New nodes
	int CreateNodes(PHCompositeNode *);

	/*
	 * fit track with SvtxTrack as input seed.
	 * \param intrack Input SvtxTrack
	 * \param invertex Input Vertex, if fit track as a primary vertex
	 */
	std::shared_ptr<PHGenFit::Track> ReFitTrack(PHCompositeNode *, const SvtxTrack* intrack, const SvtxVertex* invertex = NULL);//rcc hack: ,const bool use_svtx=true, const bool use_intt=true, const bool use_mvtx=true);
	std::shared_ptr<PHGenFit::Track> FitG4Track(PHCompositeNode *, const SvtxTrack* intrack, const SvtxVertex* invertex = NULL);//rcc hack: ,const bool use_svtx=true, const bool use_intt=true, const bool use_mvtx=true);

	//! Make SvtxTrack from PHGenFit::Track and SvtxTrack
	std::shared_ptr<SvtxTrack> MakeSvtxTrack(const SvtxTrack* svtxtrack, const std::shared_ptr<PHGenFit::Track>& genfit_track, const SvtxVertex * vertex = NULL);

	//! Fill SvtxVertexMap from GFRaveVertexes and Tracks
	bool FillSvtxVertexMap(
			const std::vector<genfit::GFRaveVertex*> & rave_vertices,
			const std::vector<genfit::Track*> & gf_tracks);

	bool pos_cov_uvn_to_rz(
			const TVector3& u,
			const TVector3& v,
			const TVector3& n,
			const TMatrixF& pos_in,
			const TMatrixF& cov_in,
			TMatrixF & pos_out,
			TMatrixF & cov_out
			) const;

	bool get_vertex_error_uvn(
			const TVector3& u,
			const TVector3& v,
			const TVector3& n,
			const TMatrixF& cov_in,
			TMatrixF & cov_out
			) const;

	bool pos_cov_XYZ_to_RZ(
			const TVector3& n,
			const TMatrixF& pos_in,
			const TMatrixF& cov_in,
			TMatrixF & pos_out,
			TMatrixF & cov_out
			) const;

	bool extrapolateTrackToRadiusPhiRZ(
			const float radius,
			std::shared_ptr<PHGenFit::Track>& rf_phgf_track,
			TVector3& pos_out,
			TMatrixF& cov_out);

	TVector3 getClusterPosAtRadius(const float radius, const SvtxTrack* intrack);
	TVector3 getClosestG4HitPos(const TVector3 target, PHCompositeNode * topNode);
	TVector3 getClosestG4HitPos(const TVector3 target, PHCompositeNode * topNode, int& nhits);
  


	/*!
	 * Get 3D Rotation Matrix that rotates frame (x,y,z) to (x',y',z')
	 * Default rotate local to global, or rotate vector in global to local representation
	 */
	TMatrixF get_rotation_matrix(
			const TVector3 x,
			const TVector3 y,
			const TVector3 z,
			const TVector3 xp = TVector3(1.,0.,0.),
			const TVector3 yp = TVector3(0.,1.,0.),
			const TVector3 zp = TVector3(0.,0.,1.)
			) const;

	//!flags
	unsigned int _flags;

	//bool _make_separate_nodes;
	OutPutMode _output_mode;

	bool _over_write_svtxtrackmap;
	bool _over_write_svtxvertexmap;

	bool _fit_primary_tracks;

	//!
	bool _use_truth_vertex;


	PHGenFit::Fitter* _fitter;

	//! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
	std::string _track_fitting_alg_name;

	int _primary_pid_guess;
	double _fit_min_pT;
	double _vertex_min_ndf;

	genfit::GFRaveVertexFactory* _vertex_finder;

	//! https://rave.hepforge.org/trac/wiki/RaveMethods
	std::string _vertexing_method;

	//PHRaveVertexFactory* _vertex_finder;

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;
	SvtxClusterMap* _clustermap;
	SvtxTrackMap* _trackmap;
	SvtxVertexMap* _vertexmap;

	//! Output Node pointers
	SvtxTrackMap* _trackmap_refit;
	SvtxTrackMap* _primary_trackmap;
	SvtxVertexMap* _vertexmap_refit;

	//! Evaluation
	//! switch eval out
	bool _do_eval;

	//! eval output filename
	std::string _eval_outname;

	TTree* _eval_tree;
	TClonesArray* _tca_particlemap;
	TClonesArray* _tca_vtxmap;
	TClonesArray* _tca_trackmap;
	TClonesArray* _tca_vertexmap;
	TClonesArray* _tca_trackmap_refit;
	TClonesArray* _tca_primtrackmap;
	TClonesArray* _tca_vertexmap_refit;

	TTree* _cluster_eval_tree;
	float _cluster_eval_tree_x;
	float _cluster_eval_tree_y;
	float _cluster_eval_tree_z;
	float _cluster_eval_tree_gx;
	float _cluster_eval_tree_gy;
	float _cluster_eval_tree_gz;

	TTree* _lost_hit_eval;
	float _lost_hit_eval_r;
	float _lost_hit_eval_x;
	float _lost_hit_eval_y;
	float _lost_hit_eval_z;
	bool _lost_hit_eval_found;
	bool _lost_hit_eval_has_svtx;
	bool _lost_hit_eval_has_intt;
	bool _lost_hit_eval_has_maps;
	bool _lost_hit_eval_in_svtx;
	bool _lost_hit_eval_in_intt;
	bool _lost_hit_eval_in_maps;
	
	TTree *_kalman_extrapolation_eval_tree;

	//data from the truth container:
	float _kalman_extrapolation_eval_tree_true_pti;
       	float _kalman_extrapolation_eval_tree_true_pxi;
	float _kalman_extrapolation_eval_tree_true_pyi;
	float _kalman_extrapolation_eval_tree_true_pzi;
	int _kalman_extrapolation_eval_tree_true_npart;
	
	//data from the svtx track alone:
	float _kalman_extrapolation_eval_tree_pti;
	float _kalman_extrapolation_eval_tree_pxi;
	float _kalman_extrapolation_eval_tree_pyi;
	float _kalman_extrapolation_eval_tree_pzi;
	int _kalman_extrapolation_eval_tree_nintt;
	int _kalman_extrapolation_eval_tree_nmvtx;

	//data from the g4 hits track:
	bool _kalman_extrapolation_eval_tree_has_g4_track;
	int _kalman_extrapolation_eval_tree_g4_nhits;

	//innermost g4tpc hit (found by looking for 30cm directly)
	int _kalman_extrapolation_eval_tree_g4_ng4hits_30_true;
	float _kalman_extrapolation_eval_tree_g4_okay_30_true;
	float _kalman_extrapolation_eval_tree_g4_phi_30_true;
	float _kalman_extrapolation_eval_tree_g4_z_30_true;
	float _kalman_extrapolation_eval_tree_g4_r_30_true;
	//g4track extrapolated to that innermost g4TPC hit.
	bool _kalman_extrapolation_eval_tree_g4_okay_ex_g4;
	float _kalman_extrapolation_eval_tree_g4_phi_ex_g4;
	float _kalman_extrapolation_eval_tree_g4_z_ex_g4;
	float _kalman_extrapolation_eval_tree_g4_r_ex_g4;
	//g4track momentum from last hit on track:
	float _kalman_extrapolation_eval_tree_g4_pt;
	float _kalman_extrapolation_eval_tree_g4_px;
	float _kalman_extrapolation_eval_tree_g4_py;
	float _kalman_extrapolation_eval_tree_g4_pz;

	//data from the cluster track:
	bool _kalman_extrapolation_eval_tree_has_cluster_track;
	int _kalman_extrapolation_eval_tree_cluster_nhits;
	//cluster track extrapolated to the ifc radius.
	bool _kalman_extrapolation_eval_tree_okay_ifc;
	float _kalman_extrapolation_eval_tree_phi_ifc;
	float _kalman_extrapolation_eval_tree_z_ifc;
	float _kalman_extrapolation_eval_tree_r_ifc;
	//with covariance info thoroughness:
	float _kalman_extrapolation_eval_tree_sigma_r_ifc;
	float _kalman_extrapolation_eval_tree_sigma_rphi_ifc;
	float _kalman_extrapolation_eval_tree_sigma_z_ifc;
	float _kalman_extrapolation_eval_tree_sigma_r_rphi_ifc;
	float _kalman_extrapolation_eval_tree_sigma_rphi_z_ifc;
	float _kalman_extrapolation_eval_tree_sigma_z_r_ifc;
	//TPC cluster position nearest R=30:
	bool _kalman_extrapolation_eval_tree_okay_30_clust;
	float _kalman_extrapolation_eval_tree_phi_30_clust;
	float _kalman_extrapolation_eval_tree_z_30_clust;
	float _kalman_extrapolation_eval_tree_r_30_clust;
	//cluster track extrapolated to the TPC cluster position.
	bool _kalman_extrapolation_eval_tree_okay_30;
	float _kalman_extrapolation_eval_tree_phi_30;
	float _kalman_extrapolation_eval_tree_z_30;
	float _kalman_extrapolation_eval_tree_r_30;
	//g4TPC hit nearest R=30:
	int _kalman_extrapolation_eval_tree_ng4hits_30_true; //number of hits nearby
	bool _kalman_extrapolation_eval_tree_okay_30_true;
	float _kalman_extrapolation_eval_tree_phi_30_true;
	float _kalman_extrapolation_eval_tree_z_30_true;
	float _kalman_extrapolation_eval_tree_r_30_true;
	//cluster extrapolated to that R=30 g4TPC hit.
	bool _kalman_extrapolation_eval_tree_okay_ex_g4;
	float _kalman_extrapolation_eval_tree_phi_ex_g4;
	float _kalman_extrapolation_eval_tree_z_ex_g4;
	float _kalman_extrapolation_eval_tree_r_ex_g4;
	//with covariance info thoroughness:
	float _kalman_extrapolation_eval_tree_sigma_r_ex_g4;
	float _kalman_extrapolation_eval_tree_sigma_rphi_ex_g4;
	float _kalman_extrapolation_eval_tree_sigma_z_ex_g4;
	float _kalman_extrapolation_eval_tree_sigma_r_rphi_ex_g4;
	float _kalman_extrapolation_eval_tree_sigma_rphi_z_ex_g4;
	float _kalman_extrapolation_eval_tree_sigma_z_r_ex_g4;
	//g4TPC hit nearest R=80:
	int _kalman_extrapolation_eval_tree_ng4hits_80_true; //number of hits nearby
	bool _kalman_extrapolation_eval_tree_okay_80_true;
	float _kalman_extrapolation_eval_tree_phi_80_true;
	float _kalman_extrapolation_eval_tree_z_80_true;
	float _kalman_extrapolation_eval_tree_r_80_true;
	//cluster extrapolated to that R=80 g4TPC hit.
	bool _kalman_extrapolation_eval_tree_okay_ex_80;
	float _kalman_extrapolation_eval_tree_phi_ex_80;
	float _kalman_extrapolation_eval_tree_z_ex_80;
	float _kalman_extrapolation_eval_tree_r_ex_80;
	//cluster track momentum from last hit on track:
	bool _kalman_extrapolation_eval_tree_p_okay;
	float _kalman_extrapolation_eval_tree_pt;
	float _kalman_extrapolation_eval_tree_px;
	float _kalman_extrapolation_eval_tree_py;
	float _kalman_extrapolation_eval_tree_pz;
	

	
	bool _do_evt_display;

};

#endif //* __PHG4TrackKalmanFitter_H__ *//

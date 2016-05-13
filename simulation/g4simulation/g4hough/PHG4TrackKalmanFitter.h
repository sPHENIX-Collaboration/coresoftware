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
	//! Default constructor
	PHG4TrackKalmanFitter(const std::string &name = "PHG4TrackKalmanFitter");

	//! dtor
	~PHG4TrackKalmanFitter();

	//!Initialization, called for initialization
	int Init(PHCompositeNode *);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *);

	//!End, write and close files
	int End(PHCompositeNode *);

	/// set verbosity
	void Verbosity(int verb) {
		verbosity = verb; // SubsysReco verbosity
	}

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

	bool is_reverse_mag_field() const {
		return _reverse_mag_field;
	}

	void set_reverse_mag_field(bool reverseMagField) {
		_reverse_mag_field = reverseMagField;
	}

	float get_mag_field_re_scaling_factor() const {
		return _mag_field_re_scaling_factor;
	}

	void set_mag_field_re_scaling_factor(float magFieldReScalingFactor) {
		_mag_field_re_scaling_factor = magFieldReScalingFactor;
	}

private:

	//! Event counter
	int _event;

	//! Get all the nodes
	int GetNodes(PHCompositeNode *);

	//!Create New nodes
	int CreateNodes(PHCompositeNode *);

	//! Refit SvtxTrack
	PHGenFit::Track* ReFitTrack(const SvtxTrack*);

	//! Make SvtxTrack from PHGenFit::Track and SvtxTrack
	SvtxTrack* MakeSvtxTrack(const SvtxTrack*, const PHGenFit::Track*);

	//! Fill SvtxVertexMap from GFRaveVertexes and Tracks
	bool FillSvtxVertexMap(
			const std::vector<genfit::GFRaveVertex*> & rave_vertices,
			const std::vector<genfit::Track*> & gf_tracks);

	//!flags
	unsigned int _flags;

	//! rescale mag field, modify the original mag field read in
	float _mag_field_re_scaling_factor;

	//! Switch to reverse Magnetic field
	bool _reverse_mag_field;

	PHGenFit::Fitter* _fitter;
	genfit::GFRaveVertexFactory* _vertex_finder;
	//PHRaveVertexFactory* _vertex_finder;

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;
	SvtxClusterMap* _clustermap;
	SvtxTrackMap* _trackmap;
	SvtxVertexMap* _vertexmap;

	//! Output Node pointers
	SvtxTrackMap* _trackmap_refit;
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
	TClonesArray* _tca_vertexmap_refit;

	bool _do_evt_display;

};

#endif //* __PHG4TrackKalmanFitter_H__ *//

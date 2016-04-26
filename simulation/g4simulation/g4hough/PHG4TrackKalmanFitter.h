/*!
 *  \file		PHG4TrackKalmanFitter.h
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4TrackKalmanFitter_H__
#define __PHG4TrackKalmanFitter_H__

#include <fun4all/SubsysReco.h>
#include <string>

class SvtxTrackMap;
class SvtxVertexMap;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;

//! \brief		Refit SvtxTracks with PHGenFit.
class PHG4TrackKalmanFitter: public SubsysReco {
public:
	//!Default constructor
	PHG4TrackKalmanFitter(const std::string &name = "PHG4TrackKalmanFitter");

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

	//Change output filename
	void set_filename(const char* file) {
		if (file)
			_eval_outname = file;
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

	//User modules
	void fill_tree(PHCompositeNode*);
	void reset_variables();

	bool is_do_eval() const {
		return _do_eval;
	}

	void set_do_eval(bool doEval) {
		_do_eval = doEval;
	}

private:
	//! switch eval out
	bool _do_eval;

	//!output filename
	std::string _eval_outname;

	//Event counter
	int _event;

	//!Get all the nodes
	int GetNodes(PHCompositeNode *);

	//!Create New nodes
	int CreateNodes(PHCompositeNode *topNode);

	//!flags
	unsigned int _flags;

	//TTrees
	TTree* _eval_tree;
	int event;
	//-- truth
	int gtrackID;
	int gflavor;
	float gpx;
	float gpy;
	float gpz;
	float gvx;
	float gvy;
	float gvz;
	//-- reco
	int trackID;
	int charge;
	int nhits;
	float px;
	float py;
	float pz;
	float dca2d;
	//-- clusters
	int clusterID[7];
	int layer[7];
	float x[7];
	float y[7];
	float z[7];
	float size_dphi[7];
	float size_dz[7];

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;
	SvtxClusterMap* _clustermap;
	SvtxTrackMap* _trackmap;

	//! Output Node pointers
	SvtxTrackMap* _trackmap_refit;
	SvtxVertexMap* _vertexmap_refit;
};

#endif //* __PHG4TrackKalmanFitter_H__ *//

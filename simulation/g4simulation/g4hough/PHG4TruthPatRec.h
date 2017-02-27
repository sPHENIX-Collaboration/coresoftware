/*!
 *  \file		PHG4TruthPatRec.h
 *  \brief		Truth Pattern Recognition
 *  \details	Using generate, GEANT level information to mimic ideal pattern recognition
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4TruthPatRec_H__
#define __PHG4TruthPatRec_H__

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>


class SvtxTrack;

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;


//! \brief		Refit SvtxTracks with PHGenFit.
class PHG4TruthPatRec: public SubsysReco {
public:

	//! Default constructor
	PHG4TruthPatRec(const std::string &name = "PHG4TruthPatRec");

	//! dtor
	~PHG4TruthPatRec();

	//!Initialization, called for initialization
	int Init(PHCompositeNode *);

	//!Initialization Run, called for initialization of a run
	int InitRun(PHCompositeNode *);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *);

	//!End, write and close files
	int End(PHCompositeNode *);

	/// set verbosity
	void Verbosity(int verb) {
		verbosity = verb; // SubsysReco verbosity
	}


private:

	//! Event counter
	int _event;

	//! Get all the nodes
	int GetNodes(PHCompositeNode *);

	//!Create New nodes
	int CreateNodes(PHCompositeNode *);

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;
	SvtxClusterMap* _clustermap;

	//! Output Node poniters
	SvtxTrackMap* _trackmap;
};

#endif //* __PHG4TruthPatRec_H__ *//

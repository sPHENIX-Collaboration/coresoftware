/*!
 *  \file		PHG4TruthPatRec.h
 *  \brief		Truth Pattern Recognition
 *  \details	Using generate, GEANT level information to mimic ideal pattern recognition
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHG4TruthPatRec_H__
#define PHG4TruthPatRec_H__

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxClusterMap;


//! \brief		Truth Pattern Recognition.
class PHG4TruthPatRec: public SubsysReco {
public:

//	enum DetectorType {MIE, MAPS_TPC, MAPS_IT_TPC, LADDER_MAPS_TPC, LADDER_MAPS_IT_TPC, LADDER_MAPS_LADDER_IT_TPC, MAPS_LADDER_IT_TPC};

	//! Default constructor
	PHG4TruthPatRec(const std::string &name = "PHG4TruthPatRec");

	//! dtor
	~PHG4TruthPatRec(){}

	//!Initialization Run, called for initialization of a run
	int InitRun(PHCompositeNode *topNode);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *topNode);

	//!End, write and close files
	int End(PHCompositeNode *topNode);

//	bool is_use_ladder_intt() const {
//		return _use_ladder_intt;
//	}
//
//	void set_use_ladder_intt(bool useLadderIntt) {
//		_use_ladder_intt = useLadderIntt;
//	}
//
//	bool is_use_ladder_maps() const {
//		return _use_ladder_maps;
//	}
//
//	void set_use_ladder_maps(bool useLadderMaps) {
//		_use_ladder_maps = useLadderMaps;
//	}

//
//	DetectorType get_detector_type() const {
//		return _detector_type;
//	}
//
//	void set_detector_type(DetectorType detectorType) {
//		_detector_type = detectorType;
//	}

private:

	//! Event counter
	int _event;

	//!Detector Type
//	DetectorType _detector_type;

//	bool _use_ladder_maps;
//	bool _use_ladder_intt;


	//! Get all the nodes
	int GetNodes(PHCompositeNode *topNode);

	//!Create New nodes
	int CreateNodes(PHCompositeNode *topNode);

	unsigned int _min_clusters_per_track;

	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;
	SvtxClusterMap* _clustermap;

	//! Output Node poniters
	SvtxTrackMap* _trackmap;
};

#endif //* __PHG4TruthPatRec_H__ *//

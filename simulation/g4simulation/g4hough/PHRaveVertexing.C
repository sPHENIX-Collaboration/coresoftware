/*!
 *  \file		PHRaveVertexing.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHRaveVertexing.h"
#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxTrackState_v1.h"
#include "SvtxHit_v1.h"
#include "SvtxHitMap.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxVertexMap_v1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>
#include <phgenfit/SpacepointMeasurement.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>

#include <GenFit/FieldManager.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <GenFit/KalmanFitterInfo.h>

//Rave
#include <rave/Version.h>
#include <rave/Track.h>
#include <rave/VertexFactory.h>
#include <rave/ConstantMagneticField.h>

//GenFit
#include <GenFit/GFRaveConverters.h>

#include <TClonesArray.h>
#include <TMatrixDSym.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TRotation.h>



#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <memory>
#include <algorithm>


#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl

#define WILD_FLOAT -9999

//#define _DEBUG_

using namespace std;

//namespace genfit {
//
//	class PHRaveVertexFactory : public GFRaveVertexFactory {
//		 void findVertices ( std::vector <  genfit::GFRaveVertex* > *, const std::vector < SvtxTrack* > &, bool use_beamspot=false ) {
//
//			  clearMap();
//
//			  try{
//			    RaveToGFVertices(GFvertices,
//			                     factory_->create(GFTracksToTracks(GFTracks, &GFStates, IdGFTrackStateMap_, 0),
//			                                      use_beamspot),
//			                     IdGFTrackStateMap_);
//			  }
//			  catch(Exception & e){
//			    clearMap();
//			    std::cerr << e.what();
//			  }
//
//			  clearMap();
//		 }
//	};
//
//}


//
//
//class PHRaveVertexFactory {
//
//public:
//	//! ctor
//	PHRaveVertexFactory(const int verbosity) {
//		rave::ConstantMagneticField mfield(0., 0., 0.); // RAVE use Tesla
//		_factory = new rave::VertexFactory(mfield, rave::VacuumPropagator(),
//				"default", verbosity);
//
//		IdGFTrackStateMap_.clear();
//	}
//
//	//! dotr
//	~PHRaveVertexFactory() {
//		clearMap();
//
//		delete _factory;
//	}
//
//	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
//			const std::vector<genfit::Track*>& tracks, const bool use_beamspot =
//					false) {
//
//		clearMap();
//
//		try {
//			genfit::RaveToGFVertices(vertices,
//					_factory->create(
//							genfit::GFTracksToTracks(tracks, NULL,
//									IdGFTrackStateMap_, 0), use_beamspot),
//					IdGFTrackStateMap_);
//		} catch (genfit::Exception & e) {
//			std::cerr << e.what();
//		}
//	}
//
//	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
//			const std::vector<genfit::Track*>& tracks,
//			std::vector<genfit::MeasuredStateOnPlane*> & GFStates,
//			const bool use_beamspot = false) {
//
//		clearMap();
//
//		try {
//			genfit::RaveToGFVertices(vertices,
//					_factory->create(
//							genfit::GFTracksToTracks(tracks, &GFStates,
//									IdGFTrackStateMap_, 0), use_beamspot),
//					IdGFTrackStateMap_);
//		} catch (genfit::Exception & e) {
//			std::cerr << e.what();
//		}
//	}
//
//private:
//	void clearMap() {
//
//		for (unsigned int i = 0; i < IdGFTrackStateMap_.size(); ++i)
//			delete IdGFTrackStateMap_[i].state_;
//
//		IdGFTrackStateMap_.clear();
//	}
//
//	std::map<int, genfit::trackAndState> IdGFTrackStateMap_;
//
//	rave::VertexFactory* _factory;
//
//};
//

/*
 * Constructor
 */
PHRaveVertexing::PHRaveVertexing(const string &name) :
		SubsysReco(name),
		_over_write_svtxtrackmap(true),
		_over_write_svtxvertexmap(true),
		_fitter( NULL),
		_primary_pid_guess(211),
		_vertex_min_ndf(20),
		_vertex_finder(NULL),
		_vertexing_method("avf-smoothing:1"),
		_truth_container(NULL),
		_clustermap(NULL),
		_trackmap(NULL),
		_vertexmap(NULL),
		_trackmap_refit(NULL),
		_primary_trackmap(NULL),
		_vertexmap_refit(NULL) {

	Verbosity(0);

	_event = 0;
}

/*
 * Init
 */
int PHRaveVertexing::Init(PHCompositeNode *topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int PHRaveVertexing::InitRun(PHCompositeNode *topNode) {

	CreateNodes(topNode);

	TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
	PHField * field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
	    field, "DafRef",
			"RKTrackRep", false);
	_fitter->set_verbosity(verbosity);

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	//LogDebug(genfit::FieldManager::getInstance()->getFieldVal(TVector3(0, 0, 0)).Z());

	_vertex_finder = new genfit::GFRaveVertexFactory(verbosity);
	_vertex_finder->setMethod(_vertexing_method.data());
	//_vertex_finder->setBeamspot();

	//_vertex_finder = new PHRaveVertexFactory(verbosity);

	if (!_vertex_finder) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}
/*!
 * process_event():
 *  Call user instructions for every event.
 *  This function contains the analysis structure.
 *
 */
int PHRaveVertexing::process_event(PHCompositeNode *topNode) {
	_event++;

	if(verbosity > 1)
		std::cout << PHWHERE << "Events processed: " << _event << std::endl;

	GetNodes(topNode);

	//! stands for Refit_GenFit_Tracks
	GenFitTrackMap gf_track_map;

	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		SvtxTrack* svtx_track = iter->second;
		if (!svtx_track)
			continue;

		if (!(svtx_track->get_ndf() >= _vertex_min_ndf))
			continue;

		genfit::Track * genfit_track = TranslateSvtxToGenFitTrack(svtx_track);
		if (!genfit_track)
			continue;

		gf_track_map.insert(pair<unsigned int, genfit::Track*> (iter->first, genfit_track));
	}

	//! find vertex using tracks
	std::vector<genfit::GFRaveVertex*> rave_vertices;

	vector<genfit::Track*> gf_tracks;

	for(auto iter : gf_track_map) {
		gf_tracks.push_back(iter.second);
	}

	if (gf_tracks.size() >= 2) {
		try {
			_vertex_finder->findVertices(&rave_vertices, gf_tracks);
		} catch (...) {
			if(verbosity > 1)
				std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
		}
	}

	FillSvtxVertexMap(rave_vertices, gf_track_map);

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int PHRaveVertexing::End(PHCompositeNode *topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * dtor
 */
PHRaveVertexing::~PHRaveVertexing() {
	delete _fitter;
	delete _vertex_finder;
}

int PHRaveVertexing::CreateNodes(PHCompositeNode *topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
  PHNodeIterator iter_dst(dstNode);

	// Create the SVTX node
	PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst(
			"PHCompositeNode", "SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (verbosity > 0)
			cout << "SVTX node added" << endl;
	}

	if (!(_over_write_svtxtrackmap)) {
		_trackmap_refit = new SvtxTrackMap_v1;
		PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
				_trackmap_refit, "SvtxTrackMapRefit", "PHObject");
		tb_node->addNode(tracks_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxTrackMapRefit node added" << endl;
	}

	if (_fit_primary_tracks) {
		_primary_trackmap = new SvtxTrackMap_v1;
		PHIODataNode<PHObject>* primary_tracks_node =
				new PHIODataNode<PHObject>(_primary_trackmap, "PrimaryTrackMap",
						"PHObject");
		tb_node->addNode(primary_tracks_node);
		if (verbosity > 0)
			cout << "Svtx/PrimaryTrackMap node added" << endl;
	}

	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
		tb_node->addNode(vertexes_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxVertexMapRefit node added" << endl;
	} else if (!findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap")) {
		_vertexmap = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap, "SvtxVertexMap", "PHObject");
		tb_node->addNode(vertexes_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxVertexMap node added" << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int PHRaveVertexing::GetNodes(PHCompositeNode * topNode) {
	//DST objects
	//Truth container
	_truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
			"G4TruthInfo");
	if (!_truth_container && _event < 2) {
		cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Clusters
	_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_clustermap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Tracks
	_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_trackmap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Vertices
	_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_vertexmap && _event < 2) {
		cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Output Svtx Tracks
	if (!(_over_write_svtxtrackmap)) {
		_trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
				"SvtxTrackMapRefit");
		if (!_trackmap_refit && _event < 2) {
			cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	// Output Primary Svtx Tracks
	if (_fit_primary_tracks) {
		_primary_trackmap = findNode::getClass<SvtxTrackMap>(topNode,
				"PrimaryTrackMap");
		if (!_primary_trackmap && _event < 2) {
			cout << PHWHERE << " PrimaryTrackMap node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	// Output Svtx Vertices
	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
				"SvtxVertexMapRefit");
		if (!_vertexmap_refit && _event < 2) {
			cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool PHRaveVertexing::FillSvtxVertexMap(
		const std::vector<genfit::GFRaveVertex*>& rave_vertices,
		const GenFitTrackMap & gf_track_map) {

	if(_over_write_svtxvertexmap){
		_vertexmap->clear();
	}


//	for (unsigned int ivtx = 0; ivtx < rave_vertices.size(); ++ivtx) {
//		genfit::GFRaveVertex* rave_vtx = rave_vertices[ivtx];

	for(genfit::GFRaveVertex* rave_vtx : rave_vertices) {

		if (!rave_vtx) {
			cerr << PHWHERE << endl;
			return false;
		}

		std::shared_ptr<SvtxVertex> svtx_vtx(new SvtxVertex_v1());

		svtx_vtx->set_chisq(rave_vtx->getChi2());
		svtx_vtx->set_ndof(rave_vtx->getNdf());
		svtx_vtx->set_position(0, rave_vtx->getPos().X());
		svtx_vtx->set_position(1, rave_vtx->getPos().Y());
		svtx_vtx->set_position(2, rave_vtx->getPos().Z());

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				svtx_vtx->set_error(i, j, rave_vtx->getCov()[i][j]);

		for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++) {
			//TODO improve speed
			const genfit::Track* rave_track =
					rave_vtx->getParameters(i)->getTrack();
			for(auto iter : gf_track_map) {
				if (iter.second == rave_track)
					svtx_vtx->insert_track(iter.first);
			}
		}

		if (_over_write_svtxvertexmap) {
			if (_vertexmap) {
				_vertexmap->insert(svtx_vtx.get());
			} else {
				LogError("!_vertexmap");
			}
		} else {
			if (_vertexmap_refit) {
				_vertexmap_refit->insert(svtx_vtx.get());
			} else {
				LogError("!_vertexmap_refit");
			}
		}

#ifdef _DEBUG_
		cout<<__LINE__<<endl;
		svtx_vtx->identify();
#endif

	} //loop over RAVE vertices

	return true;
}

genfit::Track* PHRaveVertexing::TranslateSvtxToGenFitTrack(SvtxTrack* svtx_track) {

	try {
		// The first state is extracted to PCA, second one is the one with measurement
		SvtxTrackState* svtx_state = (++(svtx_track->begin_states()))->second;

		TVector3 pos(svtx_state->get_x(), svtx_state->get_y(), svtx_state->get_z());
		TVector3 mom(svtx_state->get_px(), svtx_state->get_py(), svtx_state->get_pz());
		TMatrixDSym cov(6);
		for(int i=0;i<6;++i) {
			for(int j=0;j<6;++j) {
				cov[i][j] = svtx_state->get_error(i, j);
			}
		}

		genfit::AbsTrackRep * rep = new genfit::RKTrackRep(_primary_pid_guess);
		genfit::Track* genfit_track = new genfit::Track(rep, TVector3(0,0,0), TVector3(0,0,0));

		genfit::FitStatus * fs = new genfit::FitStatus();
		fs->setCharge(svtx_track->get_charge());
		fs->setChi2(svtx_track->get_chisq());
		fs->setNdf(svtx_track->get_ndf());
		fs->setIsFitted(true);
		fs->setIsFitConvergedFully(true);

		genfit_track->setFitStatus(fs, rep);

		genfit::TrackPoint *tp = new genfit::TrackPoint(genfit_track);

		genfit::KalmanFitterInfo* fi = new genfit::KalmanFitterInfo(tp, rep);
		tp->setFitterInfo(fi);

		genfit::MeasuredStateOnPlane * ms = new genfit::MeasuredStateOnPlane(rep);
		ms->setPosMomCov(pos, mom, cov);
		genfit::KalmanFittedStateOnPlane * kfs = new genfit::KalmanFittedStateOnPlane(*ms, 1., 1.);

		//TODO need to test
		fi->setForwardUpdate(kfs);

		genfit_track->insertPoint(tp);

		return genfit_track;
	} catch (...) {
		LogDebug("TranslateSvtxToGenFitTrack failed!");
	}

	return nullptr;
}

























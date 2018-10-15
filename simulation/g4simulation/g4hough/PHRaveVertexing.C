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

#include <phool/PHTimer.h>
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

//#define _DEBUG_

using namespace std;

/*
 * Constructor
 */
PHRaveVertexing::PHRaveVertexing(const string &name) :
		SubsysReco(name),
		_over_write_svtxvertexmap(false),
		_svtxvertexmaprefit_node_name("SvtxVertexMapRefit"),
		_fitter( NULL),
		_primary_pid_guess(211),
		_vertex_min_ndf(20),
		_vertex_finder(NULL),
		_vertexing_method("avf-smoothing:1"),
		_truth_container(NULL),
		_trackmap(NULL),
		_vertexmap(NULL),
		_t_translate(nullptr),
		_t_rave(nullptr)
		{

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
	_fitter->set_verbosity(Verbosity());

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	_vertex_finder = new genfit::GFRaveVertexFactory(Verbosity());
	_vertex_finder->setMethod(_vertexing_method.data());
	//_vertex_finder->setBeamspot();

	if (!_vertex_finder) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	_t_translate = new PHTimer("_t_translate");
	_t_translate->stop();

	_t_rave = new PHTimer("_t_rave");
	_t_rave->stop();

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

	if(Verbosity() > 1)
		std::cout << PHWHERE << "Events processed: " << _event << std::endl;

	GetNodes(topNode);

	//! stands for Refit_GenFit_Tracks
	GenFitTrackMap gf_track_map;
	vector<genfit::Track*> gf_tracks;
	if(Verbosity() > 1) _t_translate->restart();
	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		SvtxTrack* svtx_track = iter->second;
		if (!svtx_track)
			continue;

		if (!(svtx_track->get_ndf() >= _vertex_min_ndf))
			continue;

		//auto genfit_track = shared_ptr<genfit::Track> (TranslateSvtxToGenFitTrack(svtx_track));
		auto genfit_track = TranslateSvtxToGenFitTrack(svtx_track);
		if (!genfit_track)
			continue;
		gf_track_map.insert({genfit_track, iter->first});
		gf_tracks.push_back(const_cast<genfit::Track*> (genfit_track));
	}
	if(Verbosity() > 1) _t_translate->stop();

	if(Verbosity() > 1) _t_rave->restart();
	vector<genfit::GFRaveVertex*> rave_vertices;
	if (gf_tracks.size() >= 2) {
		try {
			_vertex_finder->findVertices(&rave_vertices, gf_tracks);
		} catch (...) {
			if(Verbosity() > 1)
				std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
		}
	}
	if(Verbosity() > 1) _t_rave->stop();
	FillSvtxVertexMap(rave_vertices, gf_track_map);

	for(auto iter : gf_track_map) delete iter.first;

	if(Verbosity() > 1) {
		std::cout << "=============== Timers: ===============" << std::endl;
		std::cout << "Event: " << _event << std::endl;
		std::cout << "Translate:                "<<_t_translate->get_accumulated_time()/1000. << " sec" <<std::endl;
		std::cout << "RAVE:                     "<<_t_rave->get_accumulated_time()/1000. << " sec" <<std::endl;
		std::cout << "=======================================" << std::endl;
	}

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
		if (Verbosity() > 0)
			cout << "SVTX node added" << endl;
	}

	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap_refit, _svtxvertexmaprefit_node_name.c_str(), "PHObject");
		tb_node->addNode(vertexes_node);
		if (Verbosity() > 0)
			cout << "Svtx/SvtxVertexMapRefit node added" << endl;
	} else if (!findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap")) {
		_vertexmap = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap, "SvtxVertexMap", "PHObject");
		tb_node->addNode(vertexes_node);
		if (Verbosity() > 0)
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

	// Output Svtx Vertices
	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
				_svtxvertexmaprefit_node_name.c_str());
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
//			for(auto iter : gf_track_map) {
//				if (iter.second == rave_track)
//					svtx_vtx->insert_track(iter.first);
//			}
			auto iter = gf_track_map.find(rave_track);
			if(iter != gf_track_map.end())
				svtx_vtx->insert_track(iter->second);
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
		//SvtxTrackState* svtx_state = (svtx_track->begin_states())->second;

		TVector3 pos(svtx_state->get_x(), svtx_state->get_y(), svtx_state->get_z());
		TVector3 mom(svtx_state->get_px(), svtx_state->get_py(), svtx_state->get_pz());
		TMatrixDSym cov(6);
		for(int i=0;i<6;++i) {
			for(int j=0;j<6;++j) {
				cov[i][j] = svtx_state->get_error(i, j);
			}
		}

#ifdef _DEBUG_
		{
			cout << "DEBUG" << __LINE__ << endl;
			cout << "path length:      " << svtx_state->get_pathlength() << endl;
			cout << "radius:           " << pos.Perp() << endl;
			cout << "DEBUG: " << __LINE__ << endl;
			for(int i=0; i<6; ++i){
				for(int j=0; j<6; ++j){
					cout << svtx_state->get_error(i, j) << "\t";
				}
				cout << endl;
			}

			cov.Print();
		}

#endif

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
#ifdef _DEBUG_
		{
			cout << "DEBUG: " << __LINE__ << endl;
			ms->Print();
			cout << "Orig: " << __LINE__ << endl;
			cov.Print();
			cout << "Translate: " << __LINE__ << endl;
			ms->get6DCov().Print();
		}
#endif
		genfit::KalmanFittedStateOnPlane * kfs = new genfit::KalmanFittedStateOnPlane(*ms, 1., 1.);

		//< Acording to the special order of using the stored states
		fi->setForwardUpdate(kfs);

		genfit_track->insertPoint(tp);

#ifdef _DEBUG_
//		{
//			cout << "DEBUG" << __LINE__ << endl;
//			TVector3 pos, mom;
//			TMatrixDSym cov;
//			genfit_track->getFittedState().getPosMomCov(pos, mom, cov);
//			pos.Print();
//			mom.Print();
//			cov.Print();
//		}
#endif

		return genfit_track;
	} catch (...) {
		LogDebug("TranslateSvtxToGenFitTrack failed!");
	}

	return nullptr;
}

























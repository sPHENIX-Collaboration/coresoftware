/*!
 *  \file		PHG4TrackKalmanFitter.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHG4TrackKalmanFitter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxTrackState_v1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>
#include <GenFit/FieldManager.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phgeom/PHGeomUtility.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "TClonesArray.h"
#include "TMatrixDSym.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "phgenfit/Track.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxVertexMap_v1.h"

#define LogDebug(exp)		std::cout<<"DEBUG: "<<__FILE__<<": "<<__LINE__<<": "<< #exp <<" : "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

#define _DEBUG_MODE_ 0

using namespace std;

//Rave
#include <rave/Version.h>
#include <rave/Track.h>
#include <rave/VertexFactory.h>
#include <rave/ConstantMagneticField.h>

//GenFit
#include <GenFit/GFRaveConverters.h>

class PHRaveVertexFactory {

public:
	//! ctor
	PHRaveVertexFactory(const int verbosity) {
		rave::ConstantMagneticField mfield(0., 0., 0.); // RAVE use Tesla
		_factory = new rave::VertexFactory(mfield, rave::VacuumPropagator(),
				"default", verbosity);

		IdGFTrackStateMap_.clear();
	}

	//! dotr
	~PHRaveVertexFactory() {
		clearMap();

		delete _factory;
	}

	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
			const std::vector<genfit::Track*>& tracks, const bool use_beamspot =
					false) {

		clearMap();

		try {
			genfit::RaveToGFVertices(vertices,
					_factory->create(
							genfit::GFTracksToTracks(tracks, NULL,
									IdGFTrackStateMap_, 0), use_beamspot),
					IdGFTrackStateMap_);
		} catch (genfit::Exception & e) {
			std::cerr << e.what();
		}
	}

	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
			const std::vector<genfit::Track*>& tracks,
			std::vector<genfit::MeasuredStateOnPlane*> & GFStates,
			const bool use_beamspot = false) {

		clearMap();

		try {
			genfit::RaveToGFVertices(vertices,
					_factory->create(
							genfit::GFTracksToTracks(tracks, &GFStates,
									IdGFTrackStateMap_, 0), use_beamspot),
					IdGFTrackStateMap_);
		} catch (genfit::Exception & e) {
			std::cerr << e.what();
		}
	}

private:
	void clearMap() {

		for (unsigned int i = 0; i < IdGFTrackStateMap_.size(); ++i)
			delete IdGFTrackStateMap_[i].state_;

		IdGFTrackStateMap_.clear();
	}

	std::map<int, genfit::trackAndState> IdGFTrackStateMap_;

	rave::VertexFactory* _factory;

};

/*
 * Constructor
 */
PHG4TrackKalmanFitter::PHG4TrackKalmanFitter(const string &name) :
		SubsysReco(name), _flags(NONE), _output_mode(OverwriteOriginalNode), _fit_primary_tracks(
				true), _mag_field_file_name(
				"/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root"), _mag_field_re_scaling_factor(
				1.4 / 1.5), _reverse_mag_field(true), _fitter( NULL), _track_fitting_alg_name(
				"KalmanFitterRefTrack"), _primary_pid_guess(211), _cut_min_pT(0.1), _vertex_finder(
		NULL), _vertexing_method("mvf"), _truth_container(
		NULL), _clustermap(NULL), _trackmap(NULL), _vertexmap(NULL), _trackmap_refit(
		NULL), _primary_trackmap(NULL), _vertexmap_refit(NULL), _do_eval(false), _eval_outname(
				"PHG4TrackKalmanFitter_eval.root"), _eval_tree(
		NULL), _tca_particlemap(NULL), _tca_vtxmap(NULL), _tca_trackmap(NULL), _tca_vertexmap(
		NULL), _tca_trackmap_refit(NULL), _tca_primtrackmap(NULL), _tca_vertexmap_refit(
		NULL), _do_evt_display(false) {
	_event = 0;
}

/*
 * Init
 */
int PHG4TrackKalmanFitter::Init(PHCompositeNode *topNode) {
	cout << PHWHERE << " Openning file " << _eval_outname << endl;

//	CreateNodes(topNode);

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int PHG4TrackKalmanFitter::InitRun(PHCompositeNode *topNode) {

	CreateNodes(topNode);

	TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
			_mag_field_file_name.data(),
			(_reverse_mag_field) ?
					-1. * _mag_field_re_scaling_factor :
					_mag_field_re_scaling_factor, _track_fitting_alg_name,
			"RKTrackRep", _do_evt_display);

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	//LogDebug(genfit::FieldManager::getInstance()->getFieldVal(TVector3(0, 0, 0)).Z());

	_vertex_finder = new genfit::GFRaveVertexFactory(verbosity);
	//_vertex_finder->setMethod("kalman-smoothing:1"); //! kalman-smoothing:1 is the defaul method
	_vertex_finder->setMethod(_vertexing_method.data());
	//_vertex_finder->setBeamspot();

	//_vertex_finder = new PHRaveVertexFactory(verbosity);

	if (!_vertex_finder) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	if (_do_eval) {
		PHTFileServer::get().open(_eval_outname, "RECREATE");
		init_eval_tree();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}
/*
 * process_event():
 *  Call user instructions for every event.
 *  This function contains the analysis structure.
 *
 */
int PHG4TrackKalmanFitter::process_event(PHCompositeNode *topNode) {
	_event++;
#if _DEBUG_MODE_ == 1
	cout << PHWHERE << "Events processed: " << _event << endl;
#else
	if (_event % 1000 == 0)
		cout << PHWHERE << "Events processed: " << _event << endl;
#endif

	GetNodes(topNode);

	//! stands for Refit_GenFit_Tracks
	vector<genfit::Track*> rf_gf_tracks;
	//vector<genfit::MeasuredStateOnPlane*> rf_gf_states;
	rf_gf_tracks.clear();

	vector<PHGenFit::Track*> rf_phgf_tracks;
	rf_phgf_tracks.clear();

	map<unsigned int, unsigned int> svtxtrack_genfittrack_map;

	if (_trackmap_refit)
		_trackmap_refit->empty();

	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		SvtxTrack* svtx_track = iter->second;
		if (!svtx_track)
			continue;
		if (!(svtx_track->get_pt() > _cut_min_pT))
			continue;

		//! stands for Refit_PHGenFit_Track
		PHGenFit::Track* rf_phgf_track = ReFitTrack(svtx_track);
#if _DEBUG_MODE_ == 1
		//rf_phgf_track->getGenFitTrack()->Print();
#endif
		if (rf_phgf_track) {
			svtxtrack_genfittrack_map[svtx_track->get_id()] =
					rf_phgf_tracks.size();
			rf_phgf_tracks.push_back(rf_phgf_track);
			rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());
		}
	}

	//! add tracks to event display
	if (_do_evt_display)
		_fitter->getEventDisplay()->addEvent(rf_gf_tracks);

	//! find vertex using tracks
	std::vector<genfit::GFRaveVertex*> rave_vertices;
	rave_vertices.clear();
	if (rf_gf_tracks.size() >= 2) {
		//_vertex_finder->findVertices(&rave_vertices,rf_gf_tracks,rf_gf_states);
		try {
			_vertex_finder->findVertices(&rave_vertices, rf_gf_tracks);
		} catch (...) {
			std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
		}
	}

	FillSvtxVertexMap(rave_vertices, rf_gf_tracks);

	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		PHGenFit::Track* rf_phgf_track = NULL;

		if (svtxtrack_genfittrack_map.find(iter->second->get_id())
				!= svtxtrack_genfittrack_map.end()) {
			unsigned int itrack =
					svtxtrack_genfittrack_map[iter->second->get_id()];
			rf_phgf_track = rf_phgf_tracks[itrack];
		}

		if (rf_phgf_track) {

			//FIXME figure out which vertex to use.
			SvtxVertex* vertex = NULL;
			if (_vertexmap_refit->size() > 0)
				vertex = _vertexmap_refit->get(0);

			//vertex = NULL; //DEBUG
			SvtxTrack* rf_track = MakeSvtxTrack(iter->second, rf_phgf_track,
					vertex);

			rf_phgf_tracks.push_back(rf_phgf_track);
			rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());

			if (_output_mode == MakeNewNode || _output_mode == DebugMode)
				if (_trackmap_refit)
					_trackmap_refit->insert(rf_track);

			if (_output_mode == OverwriteOriginalNode
					|| _output_mode == DebugMode)
				*(dynamic_cast<SvtxTrack_v1*>(iter->second)) =
						*(dynamic_cast<SvtxTrack_v1*>(rf_track));
		} else {
			if (_output_mode == OverwriteOriginalNode)
				_trackmap->erase(iter->first);
		}
	}

	/*!
	 * Fit track as primary track, This part need to be called after FillSvtxVertexMap
	 */
	if (_fit_primary_tracks && rave_vertices.size() > 0) {
		_primary_trackmap->empty();

		//FIXME figure out which vertex to use.
		SvtxVertex* vertex = _vertexmap_refit->get(0);
		if (vertex) {
			for (SvtxTrackMap::ConstIter iter = _trackmap->begin();
					iter != _trackmap->end(); ++iter) {
				SvtxTrack* svtx_track = iter->second;
				if (!svtx_track)
					continue;
				if (!(svtx_track->get_pt() > _cut_min_pT))
					continue;
				/*!
				 * rf_phgf_track stands for Refit_PHGenFit_Track
				 */
				PHGenFit::Track* rf_phgf_track = ReFitTrack(svtx_track,
						vertex);
				if (rf_phgf_track) {
					//FIXME figure out which vertex to use.
					SvtxVertex* vertex = NULL;
					if (_vertexmap_refit->size() > 0)
						vertex = _vertexmap_refit->get(0);
					SvtxTrack* rf_track = MakeSvtxTrack(svtx_track,
							rf_phgf_track, vertex);
					_primary_trackmap->insert(rf_track);
				}
			}
		} else {
			LogError("No vertex in SvtxVertexMapRefit!");
		}
	}

	if (_do_eval) {
		fill_eval_tree(topNode);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int PHG4TrackKalmanFitter::End(PHCompositeNode *topNode) {

	if (_do_eval) {
		PHTFileServer::get().cd(_eval_outname);
		_eval_tree->Write();
	}

	if (_do_evt_display)
		_fitter->displayEvent();

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * dtor
 */
PHG4TrackKalmanFitter::~PHG4TrackKalmanFitter() {
	delete _fitter;
	delete _vertex_finder;
}

/*
 * fill_eval_tree():
 */
void PHG4TrackKalmanFitter::fill_eval_tree(PHCompositeNode *topNode) {
	//! Make sure to reset all the TTree variables before trying to set them.
	reset_eval_variables();

	int i = 0;
	for (PHG4TruthInfoContainer::ConstIterator itr =
			_truth_container->GetPrimaryParticleRange().first;
			itr != _truth_container->GetPrimaryParticleRange().second; ++itr)
		new ((*_tca_particlemap)[i++]) (PHG4Particlev2)(
				*dynamic_cast<PHG4Particlev2*>(itr->second));

	i = 0;
	for (PHG4TruthInfoContainer::ConstVtxIterator itr =
			_truth_container->GetPrimaryVtxRange().first;
			itr != _truth_container->GetPrimaryVtxRange().second; ++itr)
		new ((*_tca_vtxmap)[i++]) (PHG4VtxPointv1)(
				*dynamic_cast<PHG4VtxPointv1*>(itr->second));

	i = 0;
	for (SvtxTrackMap::ConstIter itr = _trackmap->begin();
			itr != _trackmap->end(); ++itr)
		new ((*_tca_trackmap)[i++]) (SvtxTrack_v1)(
				*dynamic_cast<SvtxTrack_v1*>(itr->second));

	i = 0;
	for (SvtxVertexMap::ConstIter itr = _vertexmap->begin();
			itr != _vertexmap->end(); ++itr)
		new ((*_tca_vertexmap)[i++]) (SvtxVertex_v1)(
				*dynamic_cast<SvtxVertex_v1*>(itr->second));

	if (_trackmap_refit) {
		i = 0;
		for (SvtxTrackMap::ConstIter itr = _trackmap_refit->begin();
				itr != _trackmap_refit->end(); ++itr)
			new ((*_tca_trackmap_refit)[i++]) (SvtxTrack_v1)(
					*dynamic_cast<SvtxTrack_v1*>(itr->second));
	}

	if (_fit_primary_tracks) {
		i = 0;
		for (SvtxTrackMap::ConstIter itr = _primary_trackmap->begin();
				itr != _primary_trackmap->end(); ++itr)
			new ((*_tca_primtrackmap)[i++]) (SvtxTrack_v1)(
					*dynamic_cast<SvtxTrack_v1*>(itr->second));
	}

	if (_vertexmap_refit) {
		i = 0;
		for (SvtxVertexMap::ConstIter itr = _vertexmap_refit->begin();
				itr != _vertexmap_refit->end(); ++itr)
			new ((*_tca_vertexmap_refit)[i++]) (SvtxVertex_v1)(
					*dynamic_cast<SvtxVertex_v1*>(itr->second));
	}

	_eval_tree->Fill();

	return;
}

/*
 * init_eval_tree
 */
void PHG4TrackKalmanFitter::init_eval_tree() {
	if (!_tca_particlemap)
		_tca_particlemap = new TClonesArray("PHG4Particlev2");
	if (!_tca_vtxmap)
		_tca_vtxmap = new TClonesArray("PHG4VtxPointv1");

	if (!_tca_trackmap)
		_tca_trackmap = new TClonesArray("SvtxTrack_v1");
	if (!_tca_vertexmap)
		_tca_vertexmap = new TClonesArray("SvtxVertex_v1");
	if (!_tca_trackmap_refit)
		_tca_trackmap_refit = new TClonesArray("SvtxTrack_v1");
	if (_fit_primary_tracks)
		if (!_tca_primtrackmap)
			_tca_primtrackmap = new TClonesArray("SvtxTrack_v1");
	if (!_tca_vertexmap_refit)
		_tca_vertexmap_refit = new TClonesArray("SvtxVertex_v1");

	//! create TTree
	_eval_tree = new TTree("T", "PHG4TrackKalmanFitter Evaluation");

	_eval_tree->Branch("PrimaryParticle", _tca_particlemap);
	_eval_tree->Branch("TruthVtx", _tca_vtxmap);

	_eval_tree->Branch("SvtxTrack", _tca_trackmap);
	_eval_tree->Branch("SvtxVertex", _tca_vertexmap);
	_eval_tree->Branch("SvtxTrackRefit", _tca_trackmap_refit);
	if (_fit_primary_tracks)
		_eval_tree->Branch("PrimSvtxTrack", _tca_primtrackmap);
	_eval_tree->Branch("SvtxVertexRefit", _tca_vertexmap_refit);

}

/*
 * reset_eval_variables():
 *  Reset all the tree variables to their default values.
 *  Needs to be called at the start of every event
 */
void PHG4TrackKalmanFitter::reset_eval_variables() {
	_tca_particlemap->Clear();
	_tca_vtxmap->Clear();

	_tca_trackmap->Clear();
	_tca_vertexmap->Clear();
	_tca_trackmap_refit->Clear();
	if (_fit_primary_tracks)
		_tca_primtrackmap->Clear();
	_tca_vertexmap_refit->Clear();
}

int PHG4TrackKalmanFitter::CreateNodes(PHCompositeNode *topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Create the SVTX node
	PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (verbosity > 0)
			cout << "SVTX node added" << endl;
	}

	if (_output_mode == MakeNewNode || _output_mode == DebugMode) {
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

	_vertexmap_refit = new SvtxVertexMap_v1;
	PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
			_vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
	tb_node->addNode(vertexes_node);
	if (verbosity > 0)
		cout << "Svtx/SvtxVertexMapRefit node added" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int PHG4TrackKalmanFitter::GetNodes(PHCompositeNode * topNode) {
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
	if (_output_mode == MakeNewNode || _output_mode == DebugMode) {
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
	_vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
			"SvtxVertexMapRefit");
	if (!_vertexmap_refit && _event < 2) {
		cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * fit track with SvtxTrack as input seed.
 * \param intrack Input SvtxTrack
 * \param invertex Input Vertex, if fit track as a primary vertex
 */
PHGenFit::Track* PHG4TrackKalmanFitter::ReFitTrack(const SvtxTrack* intrack,
		const SvtxVertex* invertex) {
	if (!intrack) {
		cerr << PHWHERE << " Input SvtxTrack is NULL!" << endl;
		return NULL;
	}

	// prepare seed from input SvtxTrack
	TVector3 seed_mom(100, 0, 0);
	TVector3 seed_pos(0, 0, 0);
	TMatrixDSym seed_cov(6);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			seed_cov[i][j] = intrack->get_error(i, j);
		}
	}

	/*!
	 * mu+:	-13
	 * mu-:	13
	 * pi+:	211
	 * pi-:	-211
	 * e-:	11
	 * e+:	-11
	 */
	//TODO Add multiple TrackRep choices.
	//int pid = 211;
	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
	PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos, seed_mom,
			seed_cov);

	// Create measurements
	std::vector<PHGenFit::Measurement*> measurements;

	/*!
	 * if fit track as a primary track
	 */

//	if(invertex and verbosity >= 2)
//	{
//		LogDebug(invertex->size_tracks());
//		LogDebug(invertex->get_chisq());
//		LogDebug(invertex->get_ndof());
//		for (unsigned int i = 0; i < 3; i++)
//			for (unsigned int j = 0; j < 3; j++)
//			{
//				LogDebug(invertex->get_error(i,j));
//			}
//
//	}
	/*!
	 *
	 */
#if _DEBUG_MODE_ == 1
	if (invertex
//			and invertex->size_tracks() == 1
	) {
		TRandom3 rand(0);
		double dxy = 0.0007;	//7 um
		double dz = 0.003;//30 um

		TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
		TMatrixDSym cov(3);

		// Use smeared position instead of reco'd one.
		double x = rand.Gaus(0, dxy);
		double y = rand.Gaus(0, dxy);
		double z = rand.Gaus(0, dz);
		pos.SetXYZ(x, y, z);

		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		cov[i][j] = 0;

		cov[0][0] = dxy * dxy;
		cov[1][1] = dxy * dxy;
		cov[2][2] =dz * dz;

		PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
				pos, cov);
		measurements.push_back(meas);
	}
#else

	//! 1000 is a arbitrary number for now
	const double vertex_chi2_over_dnf_cut = 1000;
	const double vertex_cov_element_cut = 10000; //arbitrary cut cm*cm

	if (invertex and invertex->size_tracks() > 1
			and invertex->get_chisq() / invertex->get_ndof()
					< vertex_chi2_over_dnf_cut) {
		TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
		TMatrixDSym cov(3);
		cov.Zero();
		bool is_vertex_cov_sane = true;
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++) {

				cov(i, j) = invertex->get_error(i, j);

				if (i == j) {
					if (!(invertex->get_error(i, j) > 0
							and invertex->get_error(i, j)
									< vertex_cov_element_cut))
						is_vertex_cov_sane = false;
				}
			}

		if (is_vertex_cov_sane) {
			PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
					pos, cov);
			measurements.push_back(meas);
//			if(verbosity >= 2)
//			{
//				meas->getMeasurement()->Print();
//			}
		}
	}
#endif

	for (SvtxTrack::ConstClusterIter iter = intrack->begin_clusters();
			iter != intrack->end_clusters(); ++iter) {
		unsigned int cluster_id = *iter;
		SvtxCluster* cluster = _clustermap->get(cluster_id);
		if (!cluster) {
			LogError("No cluster Found!");
			continue;
		}
		//cluster->identify(); //DEBUG

		//unsigned int l = cluster->get_layer();

		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
		TVector3 n(cluster->get_x(), cluster->get_y(), 0);

		//TODO use u, v explicitly?
		PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
				cluster->get_phi_size(), cluster->get_z_size());

		//meas->getMeasurement()->Print();// DEBUG

		measurements.push_back(meas);
	}

	//TODO unsorted measurements, should use sorted ones?
	track->addMeasurements(measurements);

	/*!
	 *  Fit the track
	 *  ret code 0 means 0 error or good status
	 */
	if (_fitter->processTrack(track, false) != 0) {
		if (verbosity >= 1)
			LogWarning("Track fitting failed");
		return NULL;
	}

	return track;
}

/*
 * Make SvtxTrack from PHGenFit::Track and SvtxTrack
 */
SvtxTrack* PHG4TrackKalmanFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
		const PHGenFit::Track* phgf_track, const SvtxVertex* vertex) {

	double chi2 = phgf_track->get_chi2();
	double ndf = phgf_track->get_ndf();

	TVector3 vertex_position(0, 0, 0);
	double dvr2 = 0;
	double dvz2 = 0;

	if (vertex) {
		vertex_position.SetXYZ(vertex->get_x(), vertex->get_y(),
				vertex->get_z());
		dvr2 = vertex->get_error(0, 0) + vertex->get_error(1, 1);
		dvz2 = vertex->get_error(2, 2);
	}

	genfit::MeasuredStateOnPlane* gf_state_beam_line_ca = phgf_track->extrapolateToLine(
			vertex_position, TVector3(0., 0., 1.));
	TVector3 mom = gf_state_beam_line_ca->getMom();
	TVector3 pos = gf_state_beam_line_ca->getPos();
	TMatrixDSym cov = gf_state_beam_line_ca->get6DCov();

	//const SvtxTrack_v1* temp_track = static_cast<const SvtxTrack_v1*> (svtx_track);
	SvtxTrack_v1* out_track = new SvtxTrack_v1(
			*static_cast<const SvtxTrack_v1*>(svtx_track));

	/*!
	 *  1/p, u'/z', v'/z', u, v
	 *  u is defined as momentum X beam line at POCA of the beam line
	 *  v is alone the beam line
	 *  so u is the dca2d direction
	 */
	double u = gf_state_beam_line_ca->getState()[3];
	double v = gf_state_beam_line_ca->getState()[4];

	double du2 = gf_state_beam_line_ca->getCov()[3][3];
	double dv2 = gf_state_beam_line_ca->getCov()[4][4];

	if(gf_state_beam_line_ca) delete gf_state_beam_line_ca;

	out_track->set_dca2d(u);
	out_track->set_dca2d_error(sqrt(du2 + dvr2));

	genfit::MeasuredStateOnPlane* gf_state_vertex_ca = phgf_track->extrapolateToPoint(
				vertex_position);
//
//	LogDebug("Extrap to Vertex:");
//	gf_state_vertex_ca->Print();

	u = gf_state_vertex_ca->getState()[3];
	v = gf_state_vertex_ca->getState()[4];

	du2 = gf_state_vertex_ca->getCov()[3][3];
	dv2 = gf_state_vertex_ca->getCov()[4][4];

	if(gf_state_vertex_ca) delete gf_state_vertex_ca;

	double dca3d = sqrt(u * u + v * v);
	double dca3d_error = sqrt(du2 + dv2 + dvr2 + dvz2);

	out_track->set_dca(dca3d);
	out_track->set_dca_error(dca3d_error);

	out_track->set_chisq(chi2);
	out_track->set_ndf(ndf);
	out_track->set_charge(
			(_reverse_mag_field) ?
					-1. * phgf_track->get_charge() : phgf_track->get_charge());

	out_track->set_px(mom.Px());
	out_track->set_py(mom.Py());
	out_track->set_pz(mom.Pz());

	out_track->set_x(pos.X());
	out_track->set_y(pos.Y());
	out_track->set_z(pos.Z());

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			out_track->set_error(i, j, cov[i][j]);
		}
	}

	for (SvtxTrack::ConstClusterIter iter = svtx_track->begin_clusters();
			iter != svtx_track->end_clusters(); ++iter) {
		unsigned int cluster_id = *iter;
		SvtxCluster* cluster = _clustermap->get(cluster_id);
		if (!cluster) {
			LogError("No cluster Found!");
			continue;
		}
		//cluster->identify(); //DEBUG

		//unsigned int l = cluster->get_layer();

		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());

		double radius = pos.Pt();

		//TODO add exception handling
		genfit::MeasuredStateOnPlane* gf_state =
				phgf_track->extrapolateToCylinder(radius, TVector3(0, 0, 0),
						TVector3(0, 0, 1), 0);

		if (!gf_state) {
			if (verbosity > 1)
				LogWarning("Exrapolation failed!");
			continue;
		}

		SvtxTrackState* state = new SvtxTrackState_v1(radius);
		state->set_x(gf_state->getPos().x());
		state->set_y(gf_state->getPos().y());
		state->set_z(gf_state->getPos().z());

		state->set_px(gf_state->getMom().x());
		state->set_py(gf_state->getMom().y());
		state->set_pz(gf_state->getMom().z());

		//gf_state->getCov().Print();

		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				out_track->set_error(i, j, gf_state->get6DCov()[i][j]);
			}
		}

		out_track->insert_state(state);

//		std::cout<<"===============\n";
//		LogDebug(radius);
//		std::cout<<"---------------\n";
//		TVector3 temp_vec(state->get_x(),state->get_y(),state->get_z());
//		LogDebug(temp_vec.Pt());
//		//state->identify();
//		std::cout<<"---------------\n";
//		//out_track->get_state(radius)->identify();
	}

	return out_track;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool PHG4TrackKalmanFitter::FillSvtxVertexMap(
		const std::vector<genfit::GFRaveVertex*>& rave_vertices,
		const std::vector<genfit::Track*>& gf_tracks) {

	for (unsigned int ivtx = 0; ivtx < rave_vertices.size(); ++ivtx) {
		genfit::GFRaveVertex* rave_vtx = rave_vertices[ivtx];

		if (!rave_vtx) {
			cerr << PHWHERE << endl;
			return false;
		}

		SvtxVertex* svtx_vtx = new SvtxVertex_v1();

		svtx_vtx->set_chisq(rave_vtx->getChi2());
		svtx_vtx->set_ndof(rave_vtx->getNdf());
		svtx_vtx->set_position(0, rave_vtx->getPos().X());
		svtx_vtx->set_position(1, rave_vtx->getPos().Y());
		svtx_vtx->set_position(2, rave_vtx->getPos().Z());
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				svtx_vtx->set_error(i, j, rave_vtx->getCov()[i][j]);

		for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++) {
			//TODO Assume id's are sync'ed between _trackmap_refit and gf_tracks, need to change?
			const genfit::Track* rave_track =
					rave_vtx->getParameters(i)->getTrack();
			for (unsigned int j = 0; j < gf_tracks.size(); j++) {
				if (rave_track == gf_tracks[j])
					svtx_vtx->insert_track(i);
			}
		}

		if (_vertexmap_refit)
			_vertexmap_refit->insert(svtx_vtx);

		if (verbosity >= 2) {
			cout << PHWHERE << endl;
			svtx_vtx->Print();
			_vertexmap_refit->Print();
		}
	}

	return true;
}


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
#include <g4main/PHG4TruthInfoContainer.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "TMatrixDSym.h"
#include "TTree.h"
#include "TVector3.h"
#include "phgenfit/Track.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxVertexMap_v1.h"


using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
PHG4TrackKalmanFitter::PHG4TrackKalmanFitter(const string &name) :
		SubsysReco(name), _do_eval(false), _flags(NONE), _fitter( NULL), _vertex_finder( NULL ), _eval_tree( NULL) {
	//initialize
	_event = 0;
	_eval_outname = "PHG4TrackKalmanFitter_eval.root";

	reset_variables();
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::Init(PHCompositeNode *topNode) {
	cout << PHWHERE << " Openning file " << _eval_outname << endl;

	CreateNodes(topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	if(!_fitter)
	{
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	_vertex_finder = new genfit::GFRaveVertexFactory(verbosity);
	_vertex_finder->setMethod("kalman-smoothing:1");

	if (_do_eval) {
		PHTFileServer::get().open(_eval_outname, "RECREATE");
		// create TTree
		_eval_tree = new TTree("tracks", "Svtx Tracks");
		_eval_tree->Branch("event", &event, "event/I");
		_eval_tree->Branch("gtrackID", &gtrackID, "gtrackID/I");
		_eval_tree->Branch("gflavor", &gflavor, "gflavor/I");
		_eval_tree->Branch("gpx", &gpx, "gpx/F");
		_eval_tree->Branch("gpy", &gpy, "gpy/F");
		_eval_tree->Branch("gpz", &gpz, "gpz/F");
		_eval_tree->Branch("gvx", &gvx, "gvx/F");
		_eval_tree->Branch("gvy", &gvy, "gvy/F");
		_eval_tree->Branch("gvz", &gvz, "gvz/F");
		_eval_tree->Branch("trackID", &trackID, "trackID/I");
		_eval_tree->Branch("charge", &charge, "charge/I");
		_eval_tree->Branch("nhits", &nhits, "nhits/I");
		_eval_tree->Branch("px", &px, "px/F");
		_eval_tree->Branch("py", &py, "py/F");
		_eval_tree->Branch("pz", &pz, "pz/F");
		_eval_tree->Branch("dca2d", &dca2d, "dca2d/F");
		_eval_tree->Branch("clusterID", &clusterID, "clusterID[nhits]/I");
		_eval_tree->Branch("layer", &layer, "layer[nhits]/I");
		_eval_tree->Branch("x", &x, "x[nhits]/F");
		_eval_tree->Branch("y", &y, "y[nhits]/F");
		_eval_tree->Branch("z", &z, "z[nhits]/F");
		_eval_tree->Branch("size_dphi", &size_dphi, "size_dphi[nhits]/F");
		_eval_tree->Branch("size_dz", &size_dz, "size_dz[nhits]/F");
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::process_event(PHCompositeNode *topNode) {
	_event++;
	if (_event % 1000 == 0)
		cout << PHWHERE << "Events processed: " << _event << endl;

	GetNodes(topNode);

	//! stands for Refit_GenFit_Tracks
	vector<genfit::Track*> rf_gf_tracks;
	rf_gf_tracks.clear();

	std::vector<genfit::GFRaveVertex*> rave_vertices;
	rave_vertices.clear();


	for(SvtxTrackMap::ConstIter iter = _trackmap->begin(); iter != _trackmap->end();++iter)
	{
		//! stands for Refit_PHGenFit_Track
		PHGenFit::Track* rf_phgf_track = ReFitTrack(iter->second);
		SvtxTrack* rf_track = MakeSvtxTrack(iter->second,rf_phgf_track);
		_trackmap_refit->insert(rf_track);

		rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());
	}

	_vertex_finder->findVertices(&rave_vertices,rf_gf_tracks);


	if (_do_eval) {
		fill_tree(topNode);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::End(PHCompositeNode *topNode) {

	delete _fitter;

	delete _vertex_finder;

	if (_do_eval) {
		PHTFileServer::get().cd(_eval_outname);
		_eval_tree->Write();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void PHG4TrackKalmanFitter::fill_tree(PHCompositeNode *topNode) {
	//! Make sure to reset all the TTree variables before trying to set them.
	reset_variables();

	_eval_tree->Fill();

	return;
}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void PHG4TrackKalmanFitter::reset_variables() {
	event = -9999;

	//-- truth
	gtrackID = -9999;
	gflavor = -9999;
	gpx = -9999;
	gpy = -9999;
	gpz = -9999;
	gvx = -9999;
	gvy = -9999;
	gvz = -9999;

	//-- reco
	trackID = -9999;
	charge = -9999;
	nhits = -9999;
	px = -9999;
	py = -9999;
	pz = -9999;
	dca2d = -9999;

	//-- clusters
	for (int i = 0; i < 7; i++) {
		clusterID[i] = -9999;
		layer[i] = -9999;
		x[i] = -9999;
		y[i] = -9999;
		z[i] = -9999;
		size_dphi[i] = -9999;
		size_dz[i] = -9999;
	}

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

	_trackmap_refit = new SvtxTrackMap_v1;
	PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
			_trackmap_refit, "SvtxTrackMapRefit", "PHObject");
	tb_node->addNode(tracks_node);
	if (verbosity > 0)
		cout << "Svtx/SvtxTrackMap node added" << endl;

	_vertexmap_refit = new SvtxVertexMap_v1;
	PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
			_vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
	tb_node->addNode(vertexes_node);
	if (verbosity > 0)
		cout << "Svtx/SvtxVertexMap node added" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
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

	//! Input Svtx Clusters
	_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_clustermap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	//! Input Svtx Tracks
	_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_trackmap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	//! Output Svtx Tracks
	_trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
			"SvtxTrackMapRefit");
	if (!_trackmap_refit && _event < 2) {
		cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	//! Output Svtx Vertices
	_vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
			"SvtxVertexMapRefit");
	if (!_vertexmap_refit && _event < 2) {
		cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

PHGenFit::Track* PHG4TrackKalmanFitter::ReFitTrack(const SvtxTrack* intrack) {
	if(!intrack){
		cerr << PHWHERE << " Input SvtxTrack is NULL!"
						<< endl;
		return NULL;
	}

	//! prepare seed from input SvtxTrack
	TVector3 seed_mom(intrack->get_px(),intrack->get_py(),intrack->get_pz());
	TVector3 seed_pos(intrack->get_x(),intrack->get_y(),intrack->get_z());
	TMatrixDSym seed_cov(6);
	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			seed_cov[i][j] = intrack->get_error(i,j);
		}
	}

	//TODO Add multiple TrackRep choices.
	int pid = -13; //mu+
	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid);
	PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos,
			seed_mom, seed_cov);

	//! Create measurements
	std::vector<PHGenFit::Measurement*> measurements;

	for (SvtxTrack::ConstClusterIter iter = intrack->begin_clusters();
			iter != intrack->end_clusters(); ++iter) {
		unsigned int cluster_id = *iter;
		SvtxCluster* cluster = _clustermap->get(cluster_id);
		//unsigned int l = cluster->get_layer();

		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
		TVector3 n(cluster->get_x(), cluster->get_y(), 0);

		PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
				cluster->get_phi_size(), cluster->get_z_size());

		measurements.push_back(meas);
	}

	//TODO unsorted measurements, should use sorted ones?
	track->addMeasurements(measurements);

	//! Fit the track
	_fitter->processTrack(track, false);

	//TODO if not convered, make some noise

	return track;
}

SvtxTrack* PHG4TrackKalmanFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
		const PHGenFit::Track* phgf_track) {

	double chi2 = phgf_track->get_chi2();
	double ndf = phgf_track->get_ndf();

	genfit::MeasuredStateOnPlane* gf_state = phgf_track->extrapolateToLine(TVector3(0.,0.,0.), TVector3(0.,0.,1.));
	TVector3 mom = gf_state->getMom();
	TVector3 pos = gf_state->getPos();
	TMatrixDSym cov = gf_state->get6DCov();


	//const SvtxTrack_v1* temp_track = static_cast<const SvtxTrack_v1*> (svtx_track);
	SvtxTrack_v1* out_track = new SvtxTrack_v1(*static_cast<const SvtxTrack_v1*> (svtx_track));

	out_track->set_chisq(chi2);
	out_track->set_ndf(ndf);

	out_track->set_px(mom.Px());
	out_track->set_py(mom.Py());
	out_track->set_pz(mom.Pz());

	out_track->set_x(pos.X());
	out_track->set_y(pos.Y());
	out_track->set_z(pos.Z());

	for(int i=0;i<6;i++)
	{
		for(int j=i;j<6;j++)
		{
			out_track->set_error(i,j,cov[i][j]);
		}
	}

	return out_track;
}











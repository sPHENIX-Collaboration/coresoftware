#include "PHTrackSeeding.h"
#include "AssocInfoContainer.h"

#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxVertexMap_v1.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrackMap_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>



using namespace std;

PHTrackSeeding::PHTrackSeeding(const std::string& name) :
	SubsysReco(name),
	_cluster_map(nullptr),
	_vertex_map(nullptr),
	_track_map(nullptr),
	_assoc_container(nullptr)
{}

int PHTrackSeeding::Init(PHCompositeNode* topNode) {

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSeeding::InitRun(PHCompositeNode* topNode) {
	return Setup(topNode);
}

int PHTrackSeeding::process_event(PHCompositeNode* topNode) {
	return Process();
}

int PHTrackSeeding::End(PHCompositeNode* topNode) {
	return End();
}

int PHTrackSeeding::Setup(PHCompositeNode* topNode) {
	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = CreateNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSeeding::CreateNodes(PHCompositeNode* topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	PHNodeIterator iter_dst(dstNode);

	// Create the SVTX node
	PHCompositeNode* tb_node =
			dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
					"SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (Verbosity() > 0)
			cout << "SVTX node added" << endl;
	}

	_track_map = new SvtxTrackMap_v1;
	PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
			_track_map, "SvtxTrackMap", "PHObject");
	tb_node->addNode(tracks_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxTrackMap node added" << endl;

	_assoc_container = new AssocInfoContainer;
	PHIODataNode<PHObject>* assoc_node = new PHIODataNode<PHObject>(
			_assoc_container, "AssocInfoContainer", "PHObject");
	tb_node->addNode(assoc_node);
	if (Verbosity() > 0)
		cout << "Svtx/AssocInfoContainer node added" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSeeding::GetNodes(PHCompositeNode* topNode) {

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	_cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_cluster_map) {
		cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	_vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_vertex_map) {
		cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_track_map) {
		cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	_assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
	if (!_assoc_container) {
		cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

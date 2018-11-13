#include "PHInitVertexing.h"

#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxVertexMap_v1.h>
#include <g4hough/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>



using namespace std;

PHInitVertexing::PHInitVertexing(const std::string& name) :
	SubsysReco(name),
	_cluster_map(nullptr),
	_vertex_map(nullptr)
{}

int PHInitVertexing::Init(PHCompositeNode* topNode) {

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::InitRun(PHCompositeNode* topNode) {
	return Setup(topNode);
}

int PHInitVertexing::process_event(PHCompositeNode* topNode) {
	return Process();
}

int PHInitVertexing::End(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::Setup(PHCompositeNode* topNode) {
	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = CreateNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::CreateNodes(PHCompositeNode* topNode) {
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

	_vertex_map = new SvtxVertexMap_v1;
	PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
			_vertex_map, "SvtxVertexMap", "PHObject");
	tb_node->addNode(vertexes_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxVertexMap node added" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::GetNodes(PHCompositeNode* topNode) {

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	_cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_cluster_map) {
		cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Pull the reconstructed track information off the node tree...
	_vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_vertex_map) {
		cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

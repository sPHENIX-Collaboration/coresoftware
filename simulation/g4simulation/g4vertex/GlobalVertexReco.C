
#include "GlobalVertexReco.h"

#include "GlobalVertexMap_v1.h"
#include "GlobalVertex_v1.h"

#include "BbcVertexMap.h"
#include "BbcVertex.h"

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/recoConsts.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <cmath>
#include <float.h>

using namespace std;

GlobalVertexReco::GlobalVertexReco(const string &name)
  : SubsysReco(name) {
  verbosity = 0;
}

GlobalVertexReco::~GlobalVertexReco() {}

int GlobalVertexReco::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexReco::InitRun(PHCompositeNode *topNode) {

  if (verbosity > 0) {
    cout << "======================= GlobalVertexReco::InitRun() =======================" << endl;
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int GlobalVertexReco::process_event(PHCompositeNode *topNode) {
  
  if (verbosity > 1) cout << "GlobalVertexReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  GlobalVertexMap *globalmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  if (!globalmap) {
    cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap" << endl;
    exit(-1);
  }

  // if both are available, assign x,y,z from SVTX, t from BBC  
  // if only BBC is available, pass z and t from BBC forward
  // if only SVTX is available, leave t unassigned
  
  SvtxVertexMap *svtxmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  BbcVertexMap *bbcmap = findNode::getClass<BbcVertexMap>(topNode,"BbcVertexMap");

  if (svtxmap) {
    for (SvtxVertexMap::ConstIter svtxiter = svtxmap->begin();
	 svtxiter != svtxmap->end();
	 ++svtxiter) {
      const SvtxVertex* svtx = &svtxiter->second;
      GlobalVertex* vertex = new GlobalVertex_v1();

      for (unsigned int i=0; i<3; ++i) {
	vertex->set_position(i,svtx->get_position(i));
	for (unsigned int j=i; j<3; ++j) {
	  vertex->set_error(i,j,svtx->get_error(i,j));
	}
      }
      vertex->set_chisq(svtx->get_chisq());
      vertex->set_ndof(svtx->get_ndof());
      vertex->insert_vtxids(GlobalVertex::SVTX,svtx->get_id());

      // loop over the bbcmap and assign the best time
      if (bbcmap) {
	const BbcVertex *bbc_best = NULL;
	float min_sigma = FLT_MAX;
	for (BbcVertexMap::ConstIter bbciter = bbcmap->begin();
	     bbciter != bbcmap->end();
	     ++bbciter) {
	  const BbcVertex* bbc = bbciter->second;
	  
	  float combined_error = sqrt(svtx->get_error(2,2)+pow(bbc->get_z_err(),2));
	  float sigma = fabs(svtx->get_z() - bbc->get_z()) / combined_error;
	  if (sigma < min_sigma) {
	    min_sigma = sigma;
	    bbc_best = bbc;
	  }
	}

	vertex->set_t(bbc_best->get_t());
	vertex->set_t_err(bbc_best->get_t_err());
	vertex->insert_vtxids(GlobalVertex::BBC,bbc_best->get_id());
      }

      globalmap->insert(vertex);
    }

  } else if (bbcmap) {

    for (BbcVertexMap::ConstIter bbciter = bbcmap->begin();
	 bbciter != bbcmap->end();
	 ++bbciter) {
      const BbcVertex* bbc = bbciter->second;      
      GlobalVertex* vertex = new GlobalVertex_v1();

      vertex->set_t(bbc->get_t());
      vertex->set_t_err(bbc->get_t_err());
      vertex->set_z(bbc->get_z());
      vertex->set_error(2,0,0.0);
      vertex->set_error(2,1,0.0);
      vertex->set_error(2,2,pow(bbc->get_z_err(),2));
      vertex->set_error(1,2,0.0);
      vertex->set_error(0,2,0.0);
      vertex->insert_vtxids(GlobalVertex::BBC,bbc->get_id());
      
      globalmap->insert(vertex);    
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexReco::CreateNodes(PHCompositeNode *topNode) {

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the GLOBAL stuff under a sub-node directory
  PHCompositeNode* globalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","GLOBAL"));
  if (!globalNode) {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  // create the GlobalVertexMap
  GlobalVertexMap *vertexes = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  if (!vertexes) {
    vertexes = new GlobalVertexMap_v1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(vertexes,"GlobalVertexMap","PHObject");
    globalNode->addNode(VertexMapNode);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

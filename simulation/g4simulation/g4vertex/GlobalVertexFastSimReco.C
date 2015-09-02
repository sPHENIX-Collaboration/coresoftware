
#include "GlobalVertexFastSimReco.h"

#include "GlobalVertexMap_v1.h"
#include "GlobalVertex_v1.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/recoConsts.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <cmath>

using namespace std;

GlobalVertexFastSimReco::GlobalVertexFastSimReco(const string &name)
  : SubsysReco(name),
    _x_smear(NAN),
    _y_smear(NAN),
    _z_smear(NAN),
    _t_smear(NAN),
    _rand(NULL) {
  verbosity = 0;
}

GlobalVertexFastSimReco::~GlobalVertexFastSimReco() {
  if (_rand) delete _rand;
}

int GlobalVertexFastSimReco::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexFastSimReco::InitRun(PHCompositeNode *topNode) {

  if (isnan(_x_smear) ||
      isnan(_y_smear) ||
      isnan(_z_smear) ||
      isnan(_t_smear)) {
    cout << PHWHERE << "::ERROR - smearing must be defined for (x,y,z,t) via set_?_smearing(float)" << endl;
    exit(-1);
  }
  
  recoConsts *rc = recoConsts::instance();
  int seed = rc->get_IntFlag("RANDOMSEED");

  _rand = new TRandom3(seed);
  
  if (verbosity > 0) {
    cout << "=================== GlobalVertexFastSimReco::InitRun() ====================" << endl;
    cout << " x smearing: " << _x_smear << " cm " << endl;
    cout << " y smearing: " << _y_smear << " cm " << endl;
    cout << " z smearing: " << _z_smear << " cm " << endl;
    cout << " t smearing: " << _t_smear << " cm " << endl;
    cout << " random seed: " << seed << endl;
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int GlobalVertexFastSimReco::process_event(PHCompositeNode *topNode) {
  
  if (verbosity > 1) cout << "GlobalVertexFastSimReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfo) {
    cout << PHWHERE << "::ERROR - cannot find G4TruthInfo" << endl;
    exit(-1);
  }
  
  GlobalVertexMap *vertexes = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  if (!vertexes) {
    cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap" << endl;
    exit(-1);
  }

  //---------------------
  // Smear Truth Vertexes (only one per crossing right now)
  //---------------------

  PHG4VtxPoint* point = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());

  GlobalVertex* vertex = new GlobalVertex_v1();
  
  vertex->set_x(  _rand->Gaus(point->get_x() , _x_smear) );
  vertex->set_y(  _rand->Gaus(point->get_y() , _y_smear) );
  vertex->set_z(  _rand->Gaus(point->get_z() , _z_smear) );

  vertex->set_error(0,0,_x_smear*_x_smear);
  vertex->set_error(0,1,0.0);
  vertex->set_error(0,2,0.0);

  vertex->set_error(1,0,0.0);
  vertex->set_error(1,1,_y_smear*_y_smear);
  vertex->set_error(1,2,0.0);

  vertex->set_error(2,0,0.0);
  vertex->set_error(2,1,0.0);  
  vertex->set_error(2,2,_z_smear*_z_smear);

  vertex->set_t( _rand->Gaus(point->get_t(), _t_smear) );
  vertex->set_t_err( _t_smear );

  vertexes->insert(vertex);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexFastSimReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexFastSimReco::CreateNodes(PHCompositeNode *topNode) {

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
  } else {
    cout << PHWHERE << "::ERROR - GlobalVertexMap pre-exists, but should not if running FastSim" << endl;
    exit(-1);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

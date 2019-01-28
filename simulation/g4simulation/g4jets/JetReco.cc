
#include "JetReco.h"

#include "JetInput.h"
#include "JetAlgo.h"
#include "JetMap.h"
#include "JetMapV1.h"
#include "Jet.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

// standard includes
#include <iostream>
#include <vector>

using namespace std;

JetReco::JetReco(const string &name)
  : SubsysReco(name),
    _inputs(),
    _algos(),
    _algonode(),
    _inputnode(),
    _outputs() 
{}

JetReco::~JetReco() {
  for (unsigned int i=0; i<_inputs.size(); ++i) delete _inputs[i];
  _inputs.clear();
  for (unsigned int i=0; i<_algos.size(); ++i) delete _algos[i];
  _algos.clear();
  _outputs.clear();
}

int JetReco::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::InitRun(PHCompositeNode *topNode) {
  
  if (Verbosity() > 0) {
    cout << "========================== JetReco::InitRun() =============================" << endl;
    cout << " Input Selections:" << endl;
    for (unsigned int i=0; i<_inputs.size(); ++i) _inputs[i]->identify();
    cout << " Algorithms:" << endl;
    for (unsigned int i=0; i<_algos.size(); ++i) _algos[i]->identify();
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int JetReco::process_event(PHCompositeNode *topNode) {
  
  if (Verbosity() > 1) cout << "JetReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  std::vector<Jet*> inputs; // owns memory
  for (unsigned int iselect = 0; iselect < _inputs.size(); ++iselect) {
    std::vector<Jet*> parts = _inputs[iselect]->get_input(topNode);
    for (unsigned int ipart = 0; ipart < parts.size(); ++ipart) {
      inputs.push_back(parts[ipart]);
      inputs.back()->set_id(inputs.size()-1); // unique ids ensured
    }
  }

  //---------------------------
  // Run the jet reconstruction
  //---------------------------
  for (unsigned int ialgo=0; ialgo < _algos.size(); ++ialgo) {
    std::vector<Jet*> jets = _algos[ialgo]->get_jets(inputs); // owns memory

    // send the output somewhere on the DST
    FillJetNode(topNode,ialgo,jets);
  }

  // clean up input vector
  for (unsigned int i=0;i<inputs.size();++i) delete inputs[i];
  inputs.clear();
  
  if (Verbosity() > 1) cout << "JetReco::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::CreateNodes(PHCompositeNode *topNode) {

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
    
  // Create the AntiKt node if required
  PHCompositeNode* AlgoNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",_algonode.c_str()));
  if (!AlgoNode) {
    AlgoNode = new PHCompositeNode(_algonode.c_str());
    dstNode->addNode(AlgoNode);
  }
    
  // Create the Input node if required
  PHCompositeNode* InputNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",_inputnode.c_str()));
  if (!InputNode) {
    InputNode = new PHCompositeNode(_inputnode.c_str());
    AlgoNode->addNode(InputNode);
  }

  for (unsigned i=0; i<_outputs.size(); ++i) {
    JetMap *jets = findNode::getClass<JetMap>(topNode,_outputs[i]);
    if (!jets) {
      jets = new JetMapV1();
      PHIODataNode<PHObject> *JetMapNode = new PHIODataNode<PHObject>(jets,_outputs[i].c_str(),"PHObject");
      InputNode->addNode(JetMapNode);
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void JetReco::FillJetNode(PHCompositeNode *topNode, int ipos, std::vector<Jet*> jets) {

  JetMap *jetmap = NULL;
  PHTypedNodeIterator<JetMap> jetmapiter(topNode);
  PHIODataNode<JetMap> *JetMapNode = jetmapiter.find(_outputs[ipos].c_str());
  if (!JetMapNode) {
    cout << PHWHERE << " ERROR: Can't find JetMap: " << _outputs[ipos] << endl;
    exit(-1);
  } else {
    jetmap = (JetMap*)JetMapNode->getData();
  }

  jetmap->set_algo(_algos[ipos]->get_algo());
  jetmap->set_par(_algos[ipos]->get_par());
  for (unsigned int i=0; i<_inputs.size(); ++i) {
    jetmap->insert_src(_inputs[i]->get_src());
  }
  
  for (unsigned int i=0; i<jets.size(); ++i) {
    jetmap->insert(jets[i]); // map takes ownership, sets unique id
  }

  return;
}

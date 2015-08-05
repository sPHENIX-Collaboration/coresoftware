
#include "JetReco.h"

#include "JetInput.h"
#include "JetAlgo.h"
#include "JetMap.h"
#include "Jet.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>


// standard includes
#include <iostream>
#include <vector>

using namespace std;

JetReco::JetReco(const string &name)
  : SubsysReco(name),
    _inputs(),
    _algos(),
    _outputs() {
  verbosity = 0;
}

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
  
  if (verbosity >= 0) {
    cout << "========================== JetReco::InitRun() =============================" << endl;
    cout << "Input Selections:" << endl;
    cout << "Algorithms and Outputs:" << endl;    
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int JetReco::process_event(PHCompositeNode *topNode) {
  
  if (verbosity > 0) cout << "JetReco::process_event -- entered" << endl;

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
    //WriteJetOuptut(topNode,_outputs[ialgo],jets);

    // clean up --- maybe?
    for (unsigned int i=0;i<jets.size();++i) delete jets[i];
    jets.clear();
  }

  // clean up input vector
  for (unsigned int i=0;i<inputs.size();++i) delete inputs[i];
  inputs.clear();
  
  if (verbosity > 0) cout << "JetReco::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::CreateNodes(PHCompositeNode *topNode) {

  for (unsigned int i=0; i<_outputs.size(); ++i) {
    // format DST/substring/substring/substring/object
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

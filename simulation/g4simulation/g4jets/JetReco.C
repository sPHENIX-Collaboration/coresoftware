
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
  : SubsysReco(name) {
  verbosity = 0;
}

int JetReco::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::InitRun(PHCompositeNode *topNode) {
  
  if (verbosity >= 0) {
    cout << "========================== JetReco::InitRun() =============================" << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::process_event(PHCompositeNode *topNode) {
  
  if (verbosity > 0) cout << "JetReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  std::vector<Jet*> inputs;
  for (unsigned int iselect = 0; iselect < _inputs.size(); ++iselect) {
    std::vector<Jet*> parts = _inputs[iselect]->get_input(topNode);
    for (unsigned int ipart = 0; ipart < parts.size(); ++ipart) {
      inputs.push_back(parts[ipart]);
    }
  }

  //---------------------------
  // Run the jet reconstruction
  //---------------------------
  for (unsigned int ialgo=0; ialgo < _algos.size(); ++ialgo) {
    std::vector<Jet*> jets = _algos[ialgo]->get_jets(inputs);
  }

  if (verbosity > 0) cout << "JetReco::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}


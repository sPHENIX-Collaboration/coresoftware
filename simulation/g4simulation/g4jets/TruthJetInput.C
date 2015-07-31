
#include "TruthJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

// standard includes
#include <iostream>
#include <vector>

using namespace std;

TruthJetInput::TruthJetInput()
  : verbosity(0) {
}

std::vector<Jet*> TruthJetInput::get_input(PHCompositeNode *topNode) {
  
  if (verbosity > 0) cout << "TruthJetInput::process_event -- entered" << endl;

  // Pull the reconstructed track information off the node tree...
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    return std::vector<Jet*>();
  }

  std::vector<Jet*> pseudojets;
  PHG4TruthInfoContainer::Map primary_map = truthinfo->GetPrimaryMap();
  for (PHG4TruthInfoContainer::ConstIterator iter = primary_map.begin(); 
       iter != primary_map.end(); 
       ++iter) {
    PHG4Particle *part = iter->second;

    // remove some particles... stable only, muons, neutrinos, tau
    // remove acceptance... _etamin,_etamax
    
    Jet *jet = new JetV1();
    jet->set_px(part->get_px());
    jet->set_py(part->get_py());
    jet->set_pz(part->get_pz());
    jet->set_e(part->get_e());
    jet->set_id(part->get_track_id());
    pseudojets.push_back(jet);
  }

  if (verbosity > 0) cout << "TruthJetInput::process_event -- exited" << endl;

  return pseudojets;
}

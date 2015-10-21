
#include "TruthJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

// PHENIX Geant4 includes
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

// standard includes
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

TruthJetInput::TruthJetInput(Jet::SRC input)
  : _verbosity(0),
    _input(input),
    _eta_min(-4.0),
    _eta_max(+4.0) {
}

void TruthJetInput::identify(std::ostream& os) {
  os << "   TruthJetInput: G4TruthInfo to Jet::PARTICLE" << endl;
}

std::vector<Jet*> TruthJetInput::get_input(PHCompositeNode *topNode) {
  
  if (_verbosity > 0) cout << "TruthJetInput::process_event -- entered" << endl;

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

    // remove some particles (muons, taus, neutrinos)...
    // 12 == nu_e
    // 13 == muons
    // 14 == nu_mu
    // 15 == taus
    // 16 == nu_tau
    if ((abs(part->get_pid()) >= 12) && (abs(part->get_pid()) <= 16)) continue;
    
    // remove acceptance... _etamin,_etamax
    if ((part->get_px() == 0.0) && (part->get_py() == 0.0)) continue; // avoid pt=0
    float eta = asinh(part->get_pz()/sqrt(pow(part->get_px(),2)+pow(part->get_py(),2)));
    if (eta < _eta_min) continue;
    if (eta > _eta_max) continue;
    
    Jet *jet = (Jet*)(new JetV1());
    jet->set_px(part->get_px());
    jet->set_py(part->get_py());
    jet->set_pz(part->get_pz());
    jet->set_e(part->get_e());
    jet->insert_comp(Jet::PARTICLE,part->get_track_id());
    pseudojets.push_back(jet);
  }

  if (_verbosity > 0) cout << "TruthJetInput::process_event -- exited" << endl;

  return pseudojets;
}

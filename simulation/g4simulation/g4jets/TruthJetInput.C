
#include "TruthJetInput.h"

#include "JetInput.h"
#include "Jet.h"

// standard includes
#include <iostream>
#include <vector>

// fastjet includes
// #include <fastjet/JetDefinition.hh>
// #include <fastjet/PseudoJet.hh>
// #include <fastjet/ClusterSequence.hh>
// #include <fastjet/SISConePlugin.hh>

using namespace std;

TruthJetInput::TruthJetInput()
  : verbosity(0) {
}

std::vector<Jet> TruthJetInput::get_input(PHCompositeNode *topNode) {
  
  if (verbosity > 0) cout << "TruthJetInput::process_event -- entered" << endl;

  /*  
  // Pull the reconstructed track information off the node tree...
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // run fast jet
  
  std::vector<fastjet::PseudoJet> particles;
  PHG4TruthInfoContainer::Map primary_map = truthinfo->GetPrimaryMap();
  for (PHG4TruthInfoContainer::ConstIterator iter = primary_map.begin(); 
       iter != primary_map.end(); 
       ++iter) {
    PHG4Particle *part = iter->second;
    fastjet::PseudoJet pseudojet (part->get_px(),
				  part->get_py(),
				  part->get_pz(),
				  part->get_e());
    pseudojet.set_user_index(part->get_track_id());
    particles.push_back(pseudojet);
  }

  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm,0.4,fastjet::Best);
  fastjet::ClusterSequence jetFinder(particles,jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();

  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet) {
    cout << "  "
	 << fastjets[ijet].perp() << ", "
	 << fastjets[ijet].eta()  << ", "
	 << fastjets[ijet].phi()  << endl;
  }
  */
  if (verbosity > 0) cout << "TruthJetInput::process_event -- exited" << endl;

  return std::vector<Jet>();
}

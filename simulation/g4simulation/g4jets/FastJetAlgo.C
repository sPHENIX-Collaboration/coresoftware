
#include "FastJetAlgo.h"

#include "JetInput.h"
#include "Jet.h"

// standard includes
#include <iostream>
#include <vector>

// fastjet includes
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>

using namespace std;

FastJetAlgo::FastJetAlgo()
  : verbosity(0) {
}

std::vector<Jet*> FastJetAlgo::get_jets(std::vector<Jet*> particles) {
  
  if (verbosity > 0) cout << "FastJetAlgo::process_event -- entered" << endl;

  // translate to fastjet
  std::vector<fastjet::PseudoJet> pseudojets;
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart) {    
    fastjet::PseudoJet pseudojet (particles[ipart]->get_px(),
				  particles[ipart]->get_py(),
				  particles[ipart]->get_pz(),
				  particles[ipart]->get_e());
    pseudojet.set_user_index(particles[ipart]->get_id());
    pseudojets.push_back(pseudojet);
  }

  // run fast jet  
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm,0.4,fastjet::Best);
  fastjet::ClusterSequence jetFinder(pseudojets,jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();

  // print out
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet) {
    cout << "  "
	 << fastjets[ijet].perp() << ", "
	 << fastjets[ijet].eta()  << ", "
	 << fastjets[ijet].phi()  << endl;
  }

  // translate into jet output...

  if (verbosity > 0) cout << "FastJetAlgo::process_event -- exited" << endl;
  
  return std::vector<Jet*>();
}

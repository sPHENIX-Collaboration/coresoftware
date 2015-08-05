
#include "FastJetAlgo.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

// standard includes
#include <iostream>
#include <vector>

// fastjet includes
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>

using namespace std;

FastJetAlgo::FastJetAlgo(Jet::ALGO algo, float par)
  : _verbosity(0),
    _algo(algo),
    _par(par) {
}

std::vector<Jet*> FastJetAlgo::get_jets(std::vector<Jet*> particles) {
  
  if (_verbosity > 0) cout << "FastJetAlgo::process_event -- entered" << endl;

  // translate to fastjet
  std::vector<fastjet::PseudoJet> pseudojets;
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart) {    
    fastjet::PseudoJet pseudojet (particles[ipart]->get_px(),
				  particles[ipart]->get_py(),
				  particles[ipart]->get_pz(),
				  particles[ipart]->get_e());
    pseudojet.set_user_index(ipart);
    pseudojets.push_back(pseudojet);
  }

  // run fast jet
  fastjet::JetDefinition *jetdef = NULL;
  if (_algo == Jet::ANTIKT)  jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm,_par,fastjet::Best);
  else if (_algo == Jet::KT) jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm,_par,fastjet::Best);
  else if (_algo == Jet::CAMBRIDGE) jetdef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,_par,fastjet::Best);
  else return std::vector<Jet*>();
  fastjet::ClusterSequence jetFinder(pseudojets,*jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;
  
  // translate into jet output...
  std::vector<Jet*> jets;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet) {

    Jet *jet = new JetV1();
    jet->set_px(fastjets[ijet].px());
    jet->set_py(fastjets[ijet].py());
    jet->set_pz(fastjets[ijet].pz());
    jet->set_e(fastjets[ijet].e());
    jet->set_id(ijet);

    // copy components into output jet
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (unsigned int icomp = 0; icomp < comps.size(); ++icomp) {
      Jet* particle = particles[comps[icomp].user_index()];

      for (Jet::Iter iter = particle->begin_comp();
	   iter != particle->end_comp();
	   ++iter) {
      	jet->insert_comp(iter->first,iter->second);
      }
    }
    
    jets.push_back(jet);
  }

  if (_verbosity > 0) cout << "FastJetAlgo::process_event -- exited" << endl;
  
  return jets;
}

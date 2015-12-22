#include "PHPy8FwdJetTrigger.h"

#include <Pythia8/Pythia.h>

// fastjet includes
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>

using namespace std;

//__________________________________________________________
PHPy8FwdJetTrigger::PHPy8FwdJetTrigger(const std::string &name):
  PHPy8GenTrigger(name) {
  
  _verbosity = 0;
  
  _theEtaLow = 1.0;
  _theEtaHigh = 4.0; 
  
  _minPt = 10.0; 
 
}

PHPy8FwdJetTrigger::~PHPy8FwdJetTrigger() {
  if (_verbosity > 0) PrintConfig();
}

bool PHPy8FwdJetTrigger::Apply(Pythia8::Pythia *pythia) {

  if (_verbosity > 2) {
    cout << "PHPy8FwdJetTrigger::Apply - pythia event size: "
  	 << pythia->event.size() << endl;
  }
  
  // Loop over all particles in the event
  std::vector<fastjet::PseudoJet> pseudojets;
  for (int i = 0; i < pythia->event.size(); ++i) {
	
      if (pythia->event[i].status() > 0) { //only stable particles
      
	// remove some particles (muons, taus, neutrinos)...
	// 12 == nu_e
	// 13 == muons
	// 14 == nu_mu
	// 15 == taus
	// 16 == nu_tau
	if ((abs(pythia->event[i].id()) >= 12) && (abs(pythia->event[i].id()) <= 16)) continue;
    
	// remove acceptance... _etamin,_etamax
	if ((pythia->event[i].px() == 0.0) && (pythia->event[i].py() == 0.0)) continue; // avoid pt=0
  	if ( (pythia->event[i].eta() < _theEtaLow ||
  	      pythia->event[i].eta() > _theEtaHigh)) continue;


	fastjet::PseudoJet pseudojet (pythia->event[i].px(),
				  pythia->event[i].py(),
				  pythia->event[i].pz(),
				  pythia->event[i].e());
	pseudojet.set_user_index(i);
	pseudojets.push_back(pseudojet);

      }

  }

  // Call FastJet

  double _R = 1.0; 
  fastjet::JetDefinition *jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm,_R,fastjet::Best);
  fastjet::ClusterSequence jetFinder(pseudojets,*jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  bool jetFound = false; 
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet) {
    if( sqrt(pow(fastjets[ijet].px(),2) + pow(fastjets[ijet].py(),2)) > _minPt){
      jetFound = true; 
      break; 
    }
  }

  return jetFound;
}
  
void PHPy8FwdJetTrigger::SetEtaHighLow(double etaHigh, double etaLow) {
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
}

void PHPy8FwdJetTrigger::SetMinPt(double minPt) {
  _minPt = minPt; 
}

void PHPy8FwdJetTrigger::PrintConfig() {
  cout << "---------------- PHPy8FwdJetTrigger::PrintConfig --------------------" << endl;

  cout << "   Particles EtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << endl; 
  cout << "   Minimum Jet pT: " << _minPt << " GeV/c" << endl; 
  cout << "-----------------------------------------------------------------------" << endl;
}

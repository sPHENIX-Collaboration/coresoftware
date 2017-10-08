#include "PHPy8JetTrigger.h"

#include <Pythia8/Pythia.h>

// fastjet includes
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>
#include <cassert>

using namespace std;

//__________________________________________________________
PHPy8JetTrigger::PHPy8JetTrigger(const std::string &name):
  PHPy8GenTrigger(name),
  _theEtaHigh(4.0),
  _theEtaLow(1.0),
  _minPt(10.0),
  _minZ(0.0),
  _R(1.0)
 {}

PHPy8JetTrigger::~PHPy8JetTrigger() {
  if (_verbosity > 0) PrintConfig();
}

bool PHPy8JetTrigger::Apply(Pythia8::Pythia *pythia) {

  if (_verbosity > 2) {
    cout << "PHPy8JetTrigger::Apply - pythia event size: "
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

  fastjet::JetDefinition *jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm,_R, fastjet::E_scheme,fastjet::Best);
  fastjet::ClusterSequence jetFinder(pseudojets,*jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  bool jetFound = false; 
  double max_pt = -1;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet) {

      const double pt =  sqrt(pow(fastjets[ijet].px(),2) + pow(fastjets[ijet].py(),2));

      if (pt > max_pt) max_pt = pt;

      if(pt > _minPt){
	
	if(_minZ>0.0){

	  // Loop over constituents, caluclate the z of the leading particle
	  
	  float leading_Z = 0.0; 

	  float jet_ptot = sqrt( pow(fastjets[ijet].px(),2) + 
				 pow(fastjets[ijet].py(),2) + 
				 pow(fastjets[ijet].pz(),2) ); 

	  vector<fastjet::PseudoJet> constituents = fastjets[ijet].constituents();
	  for (unsigned int j=0; j<constituents.size(); j++){
	    
	    float con_ptot = sqrt( pow(constituents[j].px(),2) + 
				   pow(constituents[j].py(),2) + 
				   pow(constituents[j].pz(),2) ); 

	    float ctheta = (constituents[j].px()*fastjets[ijet].px() + 
			    constituents[j].py()*fastjets[ijet].py() + 
			    constituents[j].pz()*fastjets[ijet].pz())/(con_ptot*jet_ptot);

	    float z_constit = con_ptot*ctheta/jet_ptot; 

	    if(z_constit>leading_Z) leading_Z = z_constit; 

	  }

	  if(leading_Z>_minZ){
            jetFound = true; 
            break;
	  }

	}
	else {
          jetFound = true; 
          break;
	}
      }
  }

  if (_verbosity > 2) {
    cout << "PHPy8JetTrigger::Apply - max_pt = "<<max_pt<<", and jetFound = "<<jetFound<<endl;
  }

  return jetFound;
}
  
void PHPy8JetTrigger::SetEtaHighLow(double etaHigh, double etaLow) {

  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;

  if (_theEtaHigh<_theEtaLow)
    {
      swap(_theEtaHigh, _theEtaLow);
    }

}

void PHPy8JetTrigger::SetMinJetPt(double minPt) {
  _minPt = minPt; 
}

void PHPy8JetTrigger::SetMinLeadingZ(double minZ) {
  _minZ = minZ; 
}

void PHPy8JetTrigger::SetJetR(double R) {
  _R = R; 
}

void PHPy8JetTrigger::PrintConfig() {
  cout << "---------------- PHPy8JetTrigger::PrintConfig --------------------" << endl;

  cout << "   Particles EtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << endl; 
  cout << "   Minimum Jet pT: " << _minPt << " GeV/c" << endl; 
  cout << "   Anti-kT Radius: " << _R << endl; 
  cout << "-----------------------------------------------------------------------" << endl;
}

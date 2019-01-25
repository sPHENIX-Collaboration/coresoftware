
#include "PHPy6GenTrigger.h"
#include "PHPy6JetTrigger.h"
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>

// fastjet includes
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>
#include <cassert>

using namespace std;

//__________________________________________________________
PHPy6JetTrigger::PHPy6JetTrigger(const std::string &name):
  PHPy6GenTrigger(name),
  _theEtaHigh(4.0),
  _theEtaLow(1.0),
  _minPt(10.0),
  _R(1.0)
 {}

PHPy6JetTrigger::~PHPy6JetTrigger() {
  if (_verbosity > 0) PrintConfig();
}

bool PHPy6JetTrigger::Apply(const HepMC::GenEvent* evt) {

  
  // Loop over all particles in the event
  int idx = 0; 
  std::vector<fastjet::PseudoJet> pseudojets;
  for ( HepMC::GenEvent::particle_const_iterator p 
	  = evt->particles_begin(); p != evt->particles_end(); ++p ){
	
    idx++; 
    if ( ((*p)->status()!=1) != 0) continue; 
      
    // remove some particles (muons, taus, neutrinos)...
    // 12 == nu_e
    // 13 == muons
    // 14 == nu_mu
    // 15 == taus
    // 16 == nu_tau
    if ((abs(((*p)->pdg_id())) >= 12) && (abs(((*p)->pdg_id())) <= 16)) continue;
    
    // acceptance... _etamin,_etamax
    if (((*p)->momentum().px() == 0.0) && ((*p)->momentum().py() == 0.0)) continue; // avoid pt=0
    if ( (((*p)->momentum().pseudoRapidity()) < _theEtaLow) ||
	  (((*p)->momentum().pseudoRapidity()) > _theEtaHigh)) continue;


    fastjet::PseudoJet pseudojet ((*p)->momentum().px(),
				  (*p)->momentum().py(),
				  (*p)->momentum().pz(),
				  (*p)->momentum().e());
    pseudojet.set_user_index(idx);
    pseudojets.push_back(pseudojet);

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
      jetFound = true; 
      break; 
    }
  }

  if (_verbosity > 2) {
    cout << "PHPy6JetTrigger::Apply - max_pt = "<<max_pt<<", and jetFound = "<<jetFound<<endl;
  }

  return jetFound;
}
  
void PHPy6JetTrigger::SetEtaHighLow(double etaHigh, double etaLow) {

  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;

  if (_theEtaHigh<_theEtaLow)
    {
      swap(_theEtaHigh, _theEtaLow);
    }

}

void PHPy6JetTrigger::SetMinJetPt(double minPt) {
  _minPt = minPt; 
}

void PHPy6JetTrigger::SetJetR(double R) {
  _R = R; 
}

void PHPy6JetTrigger::PrintConfig() {
  cout << "---------------- PHPy6JetTrigger::PrintConfig --------------------" << endl;

  cout << "   Particles EtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << endl; 
  cout << "   Minimum Jet pT: " << _minPt << " GeV/c" << endl; 
  cout << "   Anti-kT Radius: " << _R << endl; 
  cout << "-----------------------------------------------------------------------" << endl;
}

#include <PHPy6GenTrigger.h>
#include <PHPy6ForwardDiElectronTrig.h>
#include "phool/PHCompositeNode.h"
#include "phool/phool.h"
#include "phool/getClass.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>

#include <cstdlib>
#include <iostream>
using namespace std;

//___________________________________________________________________________
PHPy6ForwardDiElectronTrig::PHPy6ForwardDiElectronTrig(const std::string &name): PHPy6GenTrigger(name) 
{
  ntriggered_forward_dielectron = 0;
  nconsidered_forward_dielectron = 0;
}
	
bool PHPy6ForwardDiElectronTrig::Apply( const HepMC::GenEvent* evt )
{

  // increment counter
  ++nconsidered_forward_dielectron;
	
  // Print Out Trigger Information Once, for Posterity
  static int trig_info_printed = 0;
  if ( trig_info_printed==0 )
    {
      cout << "PyForwardDiElectronTrig -- triggering on: Electron, 1.0 < |pseudorapidity| < 5.0" << endl;
      trig_info_printed = 1;
    }

  // Check the HepMC particle list - 
  // final state e+/- within 1.0 < eta < 5.0 
  // momentum > 1.0 

  int n_found = 0; 
	
  for ( HepMC::GenEvent::particle_const_iterator p 
	  = evt->particles_begin(); p != evt->particles_end(); ++p ){
    if ( (abs((*p)->pdg_id()) == 11) && ((*p)->status()==1) && 
	 ((*p)->momentum().pseudoRapidity() > 1.0) && ((*p)->momentum().pseudoRapidity() < 5.0) && 
	 (sqrt(pow((*p)->momentum().px(),2) + pow((*p)->momentum().py(),2) + pow((*p)->momentum().pz(),2))>1.0) ) {
      n_found++; 
    }
  }
  
  if(n_found>=2) {
    ++ntriggered_forward_dielectron;
    return true; 
  }

  return false;
  
}

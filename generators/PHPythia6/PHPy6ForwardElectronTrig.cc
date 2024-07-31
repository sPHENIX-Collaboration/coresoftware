#include "PHPy6ForwardElectronTrig.h"
#include "PHPy6GenTrigger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>   // for GenParticle
#include <HepMC/SimpleVector.h>  // for FourVector

#include <cmath>  // for sqrt
#include <cstdlib>
#include <iostream>

//___________________________________________________________________________
PHPy6ForwardElectronTrig::PHPy6ForwardElectronTrig(const std::string& name)
  : PHPy6GenTrigger(name)
{
  ntriggered_forward_electron = 0;
  nconsidered_forward_electron = 0;

  n_em_required = 1;
  n_ep_required = 1;
  n_comb_required = 1;
  pt_required = 0.5;
  eta_low = 1.0;
  eta_high = 5.0;

  RequireElectron = false;
  RequirePositron = false;
  RequireOR = false;
  RequireAND = false;
  RequireCOMBO = true;
}

void PHPy6ForwardElectronTrig::PrintConfig()
{
  std::cout << std::endl;
  std::cout << "PHPy6ForwardElectronTrig Configuration: " << std::endl;
  std::cout << " >=" << n_ep_required << " e+ required" << std::endl;
  std::cout << " >=" << n_em_required << " e- required" << std::endl;
  std::cout << " >=" << n_comb_required << " combined required" << std::endl;
  std::cout << " Electron transverse momentum > " << pt_required << " GeV required" << std::endl;
  std::cout << " " << eta_low << " < eta < " << eta_high << std::endl;

  if (RequireElectron)
  {
    std::cout << " RequireElectron is set" << std::endl;
  }
  if (RequirePositron)
  {
    std::cout << " RequirePositron is set" << std::endl;
  }
  if (RequireOR)
  {
    std::cout << " RequireOR is set" << std::endl;
  }
  if (RequireAND)
  {
    std::cout << " RequireAND is set" << std::endl;
  }
  if (RequireCOMBO)
  {
    std::cout << " RequireCOMBINED is set" << std::endl;
  }

  std::cout << std::endl;
}

bool PHPy6ForwardElectronTrig::Apply(const HepMC::GenEvent* evt)
{
  // increment counter
  ++nconsidered_forward_electron;

  // Print Out Trigger Information Once, for Posterity
  static int trig_info_printed = 0;
  if (trig_info_printed == 0)
  {
    PrintConfig();
    trig_info_printed = 1;
  }

  // Check the HepMC particle list -

  unsigned int n_em_found = 0;
  unsigned int n_ep_found = 0;

  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
  {
    if ((abs((*p)->pdg_id()) == 11) && ((*p)->status() == 1) &&
        ((*p)->momentum().pseudoRapidity() > eta_low) && ((*p)->momentum().pseudoRapidity() < eta_high) &&
        (sqrt(pow((*p)->momentum().px(), 2) + pow((*p)->momentum().py(), 2)) > pt_required))
    {
      if (((*p)->pdg_id()) == 11)
      {
        n_em_found++;
      }
      if (((*p)->pdg_id()) == -11)
      {
        n_ep_found++;
      }
    }
  }

  if ((RequireOR && ((n_em_found >= n_em_required) || (n_ep_found >= n_ep_required))) ||
      (RequireElectron && (n_em_found >= n_em_required)) ||
      (RequirePositron && (n_ep_found >= n_ep_required)) ||
      (RequireAND && (n_em_found >= n_em_required) && (n_ep_found >= n_ep_required)) ||
      (RequireCOMBO && (n_em_found + n_ep_found) >= n_comb_required))
  {
    ++ntriggered_forward_electron;
    return true;
  }

  return false;
}

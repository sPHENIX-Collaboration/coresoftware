#include "PHPy6JetTrigger.h"
#include "PHPy6GenTrigger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>   // for GenParticle
#include <HepMC/SimpleVector.h>  // for FourVector

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>     // for sqrt
#include <cstdlib>   // for abs
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type
#include <utility>   // for swap
#include <vector>    // for vector

using namespace std;

//__________________________________________________________
PHPy6JetTrigger::PHPy6JetTrigger(const std::string &name)
  : PHPy6GenTrigger(name)
{
}

PHPy6JetTrigger::~PHPy6JetTrigger()
{
  if (Verbosity() > 0) PrintConfig();
}

bool PHPy6JetTrigger::Apply(const HepMC::GenEvent *evt)
{
  // Loop over all particles in the event
  int idx = 0;
  std::vector<fastjet::PseudoJet> pseudojets;
  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
  {
    idx++;
    if (((*p)->status() != 1) != 0) continue;

    // remove some particles (muons, taus, neutrinos)...
    // 12 == nu_e
    // 13 == muons
    // 14 == nu_mu
    // 15 == taus
    // 16 == nu_tau
    if ((abs(((*p)->pdg_id())) >= 12) && (abs(((*p)->pdg_id())) <= 16)) continue;

    // acceptance... _etamin,_etamax
    if (((*p)->momentum().px() == 0.0) && ((*p)->momentum().py() == 0.0)) continue;  // avoid pt=0
    if ((((*p)->momentum().pseudoRapidity()) < m_theEtaLow) ||
        (((*p)->momentum().pseudoRapidity()) > m_theEtaHigh)) continue;

    fastjet::PseudoJet pseudojet((*p)->momentum().px(),
                                 (*p)->momentum().py(),
                                 (*p)->momentum().pz(),
                                 (*p)->momentum().e());
    pseudojet.set_user_index(idx);
    pseudojets.push_back(pseudojet);
  }

  // Call FastJet

  fastjet::JetDefinition *jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm, m_R, fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  bool jetFound = false;
  double max_pt = -1;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    const double pt = sqrt(pow(fastjets[ijet].px(), 2) + pow(fastjets[ijet].py(), 2));

    if (pt > max_pt) max_pt = pt;

    vector<fastjet::PseudoJet> constituents = fastjets[ijet].constituents();
    const int nconst = constituents.size();

    if (pt > m_minPt && nconst >= m_nconst)
    {
      jetFound = true;
      break;
    }
  }

  if (Verbosity() > 2)
  {
    cout << "PHPy6JetTrigger::Apply - max_pt = " << max_pt << ", and jetFound = " << jetFound << endl;
  }

  return jetFound;
}

void PHPy6JetTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  m_theEtaHigh = etaHigh;
  m_theEtaLow = etaLow;

  if (m_theEtaHigh < m_theEtaLow)
  {
    swap(m_theEtaHigh, m_theEtaLow);
  }
}

void PHPy6JetTrigger::PrintConfig()
{
  cout << "---------------- PHPy6JetTrigger::PrintConfig --------------------" << endl;

  cout << "   Particles EtaCut:  " << m_theEtaLow << " < eta < " << m_theEtaHigh << endl;
  cout << "   Minimum Jet pT: " << m_minPt << " GeV/c" << endl;
  cout << "   Anti-kT Radius: " << m_R << endl;
  cout << "-----------------------------------------------------------------------" << endl;
}

#include "PHHerwig7JetTrigger.h"

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <cmath>     // for sqrt
#include <iostream>  // for operator<<, endl, basic_ostream
#include <utility>   // for swap
#include <vector>    // for vector

//__________________________________________________________
PHHerwig7JetTrigger::PHHerwig7JetTrigger(const std::string &name)
  : PHHerwig7GenTrigger(name)
  , _theEtaHigh(4.0)
  , _theEtaLow(1.0)
  , _minPt(10.0)
  , _minZ(0.0)
  , _R(1.0)
  , _nconst(0)
{
}

PHHerwig7JetTrigger::~PHHerwig7JetTrigger()
{
  if (Verbosity() > 0)
  {
    PrintConfig();
  }
}

bool PHHerwig7JetTrigger::Apply(HepMC::GenEvent *herwig)
{
  if (Verbosity() > 2)
  {
    std::cout << "PHHerwig7JetTrigger::Apply - HepMC event size: "
              << herwig->particles_size() << std::endl;
  }

  // Loop over all particles in the event
  std::vector<fastjet::PseudoJet> pseudojets;
  for (auto iter = herwig->particles_begin(); iter != herwig->particles_end(); ++iter)
  {
    const HepMC::GenParticle *particle = *iter;
    if (particle->status() != 1)  // if not a stable particle, skip
    {
      continue;
    }

    // PHPy8JetTrigger selection
    const int absPid = std::abs(particle->pdg_id());
    if (absPid >= 12 && absPid <= 16)  // skip muons, taus, and neutrinos
    {
      continue;
    }

    const auto &momentum = particle->momentum();
    if (momentum.px() == 0.0 && momentum.py() == 0.0)
    {
      continue;
    }

    if (momentum.eta() < _theEtaLow || momentum.eta() > _theEtaHigh)
    {
      continue;
    }

    fastjet::PseudoJet pseudojet(momentum.px(), momentum.py(), momentum.pz(), momentum.e());
    pseudojet.set_user_index(particle->barcode());
    pseudojets.push_back(pseudojet);
  }

  // Call FastJet (identical to PHPy8JetTrigger)

  fastjet::JetDefinition *jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm, _R, fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  bool jetFound = false;
  double max_pt = -1;
  for (auto &fastjet : fastjets)
  {
    const double pt = sqrt(fastjet.px() * fastjet.px() + fastjet.py() * fastjet.py());

    max_pt = std::max(pt, max_pt);

    std::vector<fastjet::PseudoJet> constituents = fastjet.constituents();
    int ijet_nconst = constituents.size();

    if (pt > _minPt && ijet_nconst >= _nconst)
    {
      if (_minZ > 0.0)
      {
        // Loop over constituents, calculate the z of the leading particle

        double leading_Z = 0.0;

        double jet_ptot = sqrt(fastjet.px() * fastjet.px() +
                               fastjet.py() * fastjet.py() +
                               fastjet.pz() * fastjet.pz());

        for (auto &constituent : constituents)
        {
          double con_ptot = sqrt(constituent.px() * constituent.px() +
                                 constituent.py() * constituent.py() +
                                 constituent.pz() * constituent.pz());

          double ctheta = (constituent.px() * fastjet.px() +
                           constituent.py() * fastjet.py() +
                           constituent.pz() * fastjet.pz()) /
                          (con_ptot * jet_ptot);

          double z_constit = con_ptot * ctheta / jet_ptot;

          leading_Z = std::max(z_constit, leading_Z);
        }

        if (leading_Z > _minZ)
        {
          jetFound = true;
          break;
        }
      }
      else
      {
        jetFound = true;
        break;
      }
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "PHHerwig7JetTrigger::Apply - max_pt = " << max_pt << ", and jetFound = " << jetFound << std::endl;
  }

  return jetFound;
}

void PHHerwig7JetTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;

  if (_theEtaHigh < _theEtaLow)
  {
    std::swap(_theEtaHigh, _theEtaLow);
  }
}

void PHHerwig7JetTrigger::SetMinJetPt(double minPt)
{
  _minPt = minPt;
}

void PHHerwig7JetTrigger::SetMinLeadingZ(double minZ)
{
  _minZ = minZ;
}

void PHHerwig7JetTrigger::SetJetR(double R)
{
  _R = R;
}

void PHHerwig7JetTrigger::SetMinNumConstituents(int nconst)
{
  _nconst = nconst;
}

void PHHerwig7JetTrigger::PrintConfig() const
{
  std::cout << "---------------- PHHerwig7JetTrigger::PrintConfig --------------------" << std::endl;

  std::cout << "   Particles EtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << std::endl;
  std::cout << "   Minimum Jet pT: " << _minPt << " GeV/c" << std::endl;
  std::cout << "   Anti-kT Radius: " << _R << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
}

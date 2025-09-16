#include "PHPy8JetTrigger.h"

#include <Pythia8/Event.h>  // for Event, Particle
#include <Pythia8/Pythia.h>

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
PHPy8JetTrigger::PHPy8JetTrigger(const std::string &name)
  : PHPy8GenTrigger(name)
  , _theEtaHigh(4.0)
  , _theEtaLow(1.0)
  , _minPt(10.0)
  , _minZ(0.0)
  , _R(1.0)
  , _nconst(0)
{
}

PHPy8JetTrigger::~PHPy8JetTrigger()
{
  if (Verbosity() > 0)
  {
    PrintConfig();
  }
}

bool PHPy8JetTrigger::Apply(Pythia8::Pythia *pythia)
{
  if (Verbosity() > 2)
  {
    std::cout << "PHPy8JetTrigger::Apply - pythia event size: "
              << pythia->event.size() << std::endl;
  }

  // Loop over all particles in the event
  std::vector<fastjet::PseudoJet> pseudojets;
  for (int i = 0; i < pythia->event.size(); ++i)
  {
    if (pythia->event[i].status() > 0)
    {  // only stable particles

      // remove some particles (muons, taus, neutrinos)...
      // 12 == nu_e
      // 13 == muons
      // 14 == nu_mu
      // 15 == taus
      // 16 == nu_tau
      if ((abs(pythia->event[i].id()) >= 12) && (abs(pythia->event[i].id()) <= 16))
      {
        continue;
      }

      // remove acceptance... _etamin,_etamax
      if ((pythia->event[i].px() == 0.0) && (pythia->event[i].py() == 0.0))
      {
        continue;  // avoid pt=0
      }
      if ((pythia->event[i].eta() < _theEtaLow ||
           pythia->event[i].eta() > _theEtaHigh))
      {
        continue;
      }

      fastjet::PseudoJet pseudojet(pythia->event[i].px(),
                                   pythia->event[i].py(),
                                   pythia->event[i].pz(),
                                   pythia->event[i].e());
      pseudojet.set_user_index(i);
      pseudojets.push_back(pseudojet);
    }
  }

  // Call FastJet

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
    std::cout << "PHPy8JetTrigger::Apply - max_pt = " << max_pt << ", and jetFound = " << jetFound << std::endl;
  }

  return jetFound;
}

void PHPy8JetTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;

  if (_theEtaHigh < _theEtaLow)
  {
    std::swap(_theEtaHigh, _theEtaLow);
  }
}

void PHPy8JetTrigger::SetMinJetPt(double minPt)
{
  _minPt = minPt;
}

void PHPy8JetTrigger::SetMinLeadingZ(double minZ)
{
  _minZ = minZ;
}

void PHPy8JetTrigger::SetJetR(double R)
{
  _R = R;
}

void PHPy8JetTrigger::SetMinNumConstituents(int nconst)
{
  _nconst = nconst;
}

void PHPy8JetTrigger::PrintConfig() const
{
  std::cout << "---------------- PHPy8JetTrigger::PrintConfig --------------------" << std::endl;

  std::cout << "   Particles EtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << std::endl;
  std::cout << "   Minimum Jet pT: " << _minPt << " GeV/c" << std::endl;
  std::cout << "   Anti-kT Radius: " << _R << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
}

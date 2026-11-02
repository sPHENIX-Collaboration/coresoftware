#include "PHHerwig7ParticleTrigger.h"

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/SimpleVector.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>

namespace
{
  bool inRange(double value, double low, double high)
  {
    return value >= low && value <= high;
  }

  void printCut(const std::string &name, double low, double high)
  {
    if (std::isfinite(low) || std::isfinite(high))
    {
      std::cout << "   " << name << " cut:  " << low << " < " << name << " < " << high << std::endl;
    }
  }

  // HepMC::FourVector has no rapidity accessor
  double rapidity(const HepMC::FourVector &mom)
  {
    return 0.5 * std::log((mom.e() + mom.pz()) / (mom.e() - mom.pz()));
  }

  bool hasHadronParent(const HepMC::GenVertex *prod)
  {
    for (auto pin = prod->particles_in_const_begin(); pin != prod->particles_in_const_end(); ++pin)
    {
      const int absId = std::abs((*pin)->pdg_id());
      if (absId > 100 && absId != 2212 && absId != 2112)
      {
        return true;
      }
    }
    return false;
  }

  bool hasListedParent(const HepMC::GenVertex *prod, const std::vector<int> &parents)
  {
    for (auto pin = prod->particles_in_const_begin(); pin != prod->particles_in_const_end(); ++pin)
    {
      const int absId = std::abs((*pin)->pdg_id());
      if (std::any_of(parents.begin(), parents.end(), [absId](int parent)
                      { return std::abs(parent) == absId; }))
      {
        return true;
      }
    }
    return false;
  }
}  // namespace

//__________________________________________________________
PHHerwig7ParticleTrigger::PHHerwig7ParticleTrigger(const std::string &name)
  : PHHerwig7GenTrigger(name)
{
}

PHHerwig7ParticleTrigger::~PHHerwig7ParticleTrigger()
{
  if (Verbosity() > 0)
  {
    PrintConfig();
  }
}

bool PHHerwig7ParticleTrigger::Apply(HepMC::GenEvent *herwig)
{
  if (Verbosity() > 2)
  {
    std::cout << "PHHerwig7ParticleTrigger::Apply - HepMC event size: "
              << herwig->particles_size() << std::endl;
  }

  for (auto iter = herwig->particles_begin(); iter != herwig->particles_end(); ++iter)
  {
    const HepMC::GenParticle *particle = *iter;
    if (std::find(m_Particles.begin(), m_Particles.end(), particle->pdg_id()) == m_Particles.end())
    {
      continue;
    }
    if (m_StableParticleOnly && particle->status() != 1)
    {
      continue;
    }

    const HepMC::FourVector &mom = particle->momentum();
    if (mom.px() == 0.0 && mom.py() == 0.0)  // beam-axis particles have no defined eta
    {
      continue;
    }
    const double eta = mom.eta();
    if (!inRange(mom.perp(), m_PtLow, m_PtHigh) ||
        !inRange(mom.rho(), m_PLow, m_PHigh) ||
        !inRange(mom.pz(), m_PzLow, m_PzHigh) ||
        !inRange(rapidity(mom), m_YLow, m_YHigh) ||
        !inRange(eta, m_EtaLow, m_EtaHigh) ||
        !inRange(std::abs(eta), m_AbsEtaLow, m_AbsEtaHigh))
    {
      continue;
    }

    const HepMC::GenVertex *prod = particle->production_vertex();
    if (m_PromptOnly && prod && hasHadronParent(prod))
    {
      continue;
    }
    if (!m_Parents.empty() && !(prod && hasListedParent(prod, m_Parents)))
    {
      continue;
    }

    if (Verbosity() > 5)
    {
      std::cout << GetName() << ": pass on id " << particle->pdg_id() << " status " << particle->status()
                << " pt " << mom.perp() << " pz " << mom.pz() << " p " << mom.rho()
                << " eta " << eta << " y " << rapidity(mom) << std::endl;
    }
    return true;
  }
  return false;
}

void PHHerwig7ParticleTrigger::AddParticles(const std::string &particles)
{
  AddParticles(convertToInts(particles));
}

void PHHerwig7ParticleTrigger::AddParticles(int particle)
{
  m_Particles.push_back(particle);
}

void PHHerwig7ParticleTrigger::AddParticles(const std::vector<int> &particles)
{
  m_Particles.insert(m_Particles.end(), particles.begin(), particles.end());
}

void PHHerwig7ParticleTrigger::AddParents(const std::string &parents)
{
  AddParents(convertToInts(parents));
}

void PHHerwig7ParticleTrigger::AddParents(int parent)
{
  m_Parents.push_back(parent);
}

void PHHerwig7ParticleTrigger::AddParents(const std::vector<int> &parents)
{
  m_Parents.insert(m_Parents.end(), parents.begin(), parents.end());
}

void PHHerwig7ParticleTrigger::SetPtHighLow(double ptHigh, double ptLow)
{
  m_PtHigh = std::max(ptHigh, ptLow);
  m_PtLow = std::min(ptHigh, ptLow);
}

void PHHerwig7ParticleTrigger::SetPHighLow(double pHigh, double pLow)
{
  m_PHigh = std::max(pHigh, pLow);
  m_PLow = std::min(pHigh, pLow);
}

void PHHerwig7ParticleTrigger::SetYHighLow(double YHigh, double YLow)
{
  m_YHigh = std::max(YHigh, YLow);
  m_YLow = std::min(YHigh, YLow);
}

void PHHerwig7ParticleTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  m_EtaHigh = std::max(etaHigh, etaLow);
  m_EtaLow = std::min(etaHigh, etaLow);
}

void PHHerwig7ParticleTrigger::SetAbsEtaHighLow(double etaHigh, double etaLow)
{
  m_AbsEtaHigh = std::max(etaHigh, etaLow);
  m_AbsEtaLow = std::min(etaHigh, etaLow);
}

void PHHerwig7ParticleTrigger::SetPzHighLow(double pzHigh, double pzLow)
{
  m_PzHigh = std::max(pzHigh, pzLow);
  m_PzLow = std::min(pzHigh, pzLow);
}

void PHHerwig7ParticleTrigger::PrintConfig() const
{
  std::cout << "---------------- PHHerwig7ParticleTrigger::PrintConfig --------------------" << std::endl;

  std::cout << (m_StableParticleOnly ? "Process stable particles only." : "Process both unstable and stable particles.") << std::endl;
  if (m_PromptOnly)
  {
    std::cout << "Reject particles with a hadron parent." << std::endl;
  }

  std::cout << "   Particles: ";
  for (int particle : m_Particles)
  {
    std::cout << particle << "  ";
  }
  std::cout << std::endl;

  std::cout << "   Parents: ";
  for (int parent : m_Parents)
  {
    std::cout << parent << "  ";
  }
  std::cout << std::endl;

  printCut("Y", m_YLow, m_YHigh);
  printCut("eta", m_EtaLow, m_EtaHigh);
  printCut("|eta|", m_AbsEtaLow, m_AbsEtaHigh);
  printCut("pT", m_PtLow, m_PtHigh);
  printCut("p", m_PLow, m_PHigh);
  printCut("pz", m_PzLow, m_PzHigh);
  std::cout << "-----------------------------------------------------------------------" << std::endl;
}

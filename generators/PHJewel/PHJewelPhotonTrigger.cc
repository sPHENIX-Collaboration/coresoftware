#include "PHJewelPhotonTrigger.h"

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace
{
  bool isNeutrino(int pid)
  {
    const int a = std::abs(pid);
    return (a == 12 || a == 14 || a == 16);
  }

  double dphiPi(double a, double b)
  {
    double d = std::fabs(a - b);
    while (d > M_PI)
    {
      d = std::fabs(d - 2.0 * M_PI);
    }
    return d;
  }

  // a photon is a decay photon if any immediate parent is a hadron
  bool isDecayPhoton(const HepMC::GenParticle *p)
  {
    const HepMC::GenVertex *pv = p->production_vertex();
    if (!pv)
    {
      return false;
    }
    for (HepMC::GenVertex::particles_in_const_iterator it = pv->particles_in_const_begin();
         it != pv->particles_in_const_end(); ++it)
    {
      const int a = std::abs((*it)->pdg_id());
      if (a > 100 && a != 2212 && a != 2112)
      {
        return true;
      }
      if (a == 111 || a == 221)
      {
        return true;
      }
    }
    return false;
  }
}  // namespace

PHJewelPhotonTrigger::PHJewelPhotonTrigger(const std::string &name)
  : PHJewelGenTrigger(name)
{
}

bool PHJewelPhotonTrigger::Apply(const HepMC::GenEvent *evt)
{
  m_LeadPhotonPt = 0.0;

  std::vector<const HepMC::GenParticle *> fs;
  for (HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it)
  {
    if ((*it)->status() == 1)
    {
      fs.push_back(*it);
    }
  }

  for (const HepMC::GenParticle *p : fs)
  {
    if (p->pdg_id() != 22)
    {
      continue;
    }
    const HepMC::FourVector &m = p->momentum();
    const double pt = m.perp();
    if (pt < m_PtMin || std::fabs(m.eta()) > m_EtaMax)
    {
      continue;
    }
    if (m_Prompt && isDecayPhoton(p))
    {
      continue;
    }
    if (m_Iso)
    {
      double cone = 0.0;
      for (const HepMC::GenParticle *q : fs)
      {
        if (q == p || isNeutrino(q->pdg_id()))
        {
          continue;
        }
        const HepMC::FourVector &mq = q->momentum();
        const double de = mq.eta() - m.eta();
        const double dp = dphiPi(mq.phi(), m.phi());
        if (de * de + dp * dp < m_IsoR * m_IsoR)
        {
          cone += mq.perp();
        }
      }
      if (cone >= m_IsoFrac * pt)
      {
        continue;
      }
    }
    m_LeadPhotonPt = std::max(pt, m_LeadPhotonPt);
  }

  const bool pass = m_LeadPhotonPt > 0.0;
  if (Verbosity() > 1)
  {
    std::cout << GetName() << ": leading accepted photon pT " << m_LeadPhotonPt << " pass " << pass << std::endl;
  }
  return pass;
}

void PHJewelPhotonTrigger::PrintConfig() const
{
  std::cout << GetName() << ": photon pT>=" << m_PtMin << " |eta|<=" << m_EtaMax
            << " isolation " << (m_Iso ? "on" : "off");
  if (m_Iso)
  {
    std::cout << " (R=" << m_IsoR << ", frac=" << m_IsoFrac << ")";
  }
  std::cout << " prompt " << (m_Prompt ? "on" : "off") << std::endl;
}

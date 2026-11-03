#include "PHJewelJetTrigger.h"

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

PHJewelJetTrigger::PHJewelJetTrigger(const std::string &name)
  : PHJewelGenTrigger(name)
{
}

bool PHJewelJetTrigger::Apply(const HepMC::GenEvent *evt)
{
  m_LeadJetPt = 0.0;

  std::vector<fastjet::PseudoJet> input;
  for (HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it)
  {
    const HepMC::GenParticle *p = *it;
    if (p->status() != 1)
    {
      continue;
    }
    const int apid = std::abs(p->pdg_id());
    if (apid == 12 || apid == 14 || apid == 16)
    {
      continue;
    }
    const HepMC::FourVector &m = p->momentum();
    input.emplace_back(m.px(), m.py(), m.pz(), m.e());
  }

  if (!input.empty())
  {
    const fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, m_R);
    fastjet::ClusterSequence cs(input, jetDef);
    const std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(1.0));
    for (const fastjet::PseudoJet &j : jets)
    {
      if (std::fabs(j.eta()) < m_EtaMax)
      {
        m_LeadJetPt = j.pt();
        break;
      }
    }
  }

  const bool pass = m_LeadJetPt >= m_PtMin && m_LeadJetPt >= m_WindowLo && m_LeadJetPt < m_WindowHi;
  if (Verbosity() > 1)
  {
    std::cout << GetName() << ": leading jet pT " << m_LeadJetPt << " pass " << pass << std::endl;
  }
  return pass;
}

void PHJewelJetTrigger::PrintConfig() const
{
  std::cout << GetName() << ": anti-kT R=" << m_R << " pT>=" << m_PtMin << " |eta|<" << m_EtaMax
            << " leading-jet window [" << m_WindowLo << ", " << m_WindowHi << ")" << std::endl;
}

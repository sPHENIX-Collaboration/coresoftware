#include "FastJetAlgo.h"

#include "Jet.h"
#include "Jetv1.h"

#include <phool/phool.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>  // for FunctionOfPse...
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

// SoftDrop includes
#include <fastjet/contrib/RecursiveSymmetryCutBase.hh>  // for RecursiveSymm...
#include <fastjet/contrib/SoftDrop.hh>

#include <TSystem.h>

// standard includes
#include <cmath>  // for isfinite
#include <iostream>
#include <map>      // for _Rb_tree_iterator
#include <memory>   // for allocator_traits<>::value_type
#include <string>   // for operator<<
#include <utility>  // for pair
#include <vector>

FastJetAlgo::FastJetAlgo(Jet::ALGO algo, float par, int verbosity)
  : m_Verbosity(verbosity)
  , m_AlgoFlag(algo)
  , m_Par(par)
  , m_SDFlag(false)
  , m_SDBeta(0)
  , m_SDZCut(0.1)
{
  fastjet::ClusterSequence clusseq;
  if (m_Verbosity > 0)
  {
    clusseq.print_banner();
  }
  else
  {
    std::ostringstream nullstream;
    clusseq.set_fastjet_banner_stream(&nullstream);
    clusseq.print_banner();
    clusseq.set_fastjet_banner_stream(&std::cout);
  }
}

void FastJetAlgo::identify(std::ostream& os)
{
  os << "   FastJetAlgo: ";
  if (m_AlgoFlag == Jet::ANTIKT)
  {
    os << "ANTIKT r=" << m_Par;
  }
  else if (m_AlgoFlag == Jet::KT)
  {
    os << "KT r=" << m_Par;
  }
  else if (m_AlgoFlag == Jet::CAMBRIDGE)
  {
    os << "CAMBRIDGE r=" << m_Par;
  }
  os << std::endl;
}

std::vector<Jet*> FastJetAlgo::get_jets(std::vector<Jet*> particles)
{
  if (m_Verbosity > 1) std::cout << "FastJetAlgo::process_event -- entered" << std::endl;

  // translate to fastjet
  std::vector<fastjet::PseudoJet> pseudojets;
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart)
  {
    // fastjet performs strangely with exactly (px,py,pz,E) =
    // (0,0,0,0) inputs, such as placeholder towers or those with
    // zero'd out energy after CS. this catch also in FastJetAlgoSub
    if (particles[ipart]->get_e() == 0.) continue;
    if (!std::isfinite(particles[ipart]->get_px()) ||
        !std::isfinite(particles[ipart]->get_py()) ||
        !std::isfinite(particles[ipart]->get_pz()) ||
        !std::isfinite(particles[ipart]->get_e()))
    {
      std::cout << PHWHERE << " invalid particle kinematics:"
                << " px: " << particles[ipart]->get_px()
                << " py: " << particles[ipart]->get_py()
                << " pz: " << particles[ipart]->get_pz()
                << " e: " << particles[ipart]->get_e() << std::endl;
      gSystem->Exit(1);
    }
    fastjet::PseudoJet pseudojet(particles[ipart]->get_px(),
                                 particles[ipart]->get_py(),
                                 particles[ipart]->get_pz(),
                                 particles[ipart]->get_e());
    pseudojet.set_user_index(ipart);
    pseudojets.push_back(pseudojet);
  }
  // run fast jet
  fastjet::JetDefinition* jetdef = nullptr;
  if (m_AlgoFlag == Jet::ANTIKT)
  {
    jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm, m_Par, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_AlgoFlag == Jet::KT)
  {
    jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, m_Par, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_AlgoFlag == Jet::CAMBRIDGE)
  {
    jetdef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, m_Par, fastjet::E_scheme, fastjet::Best);
  }
  else
  {
    return std::vector<Jet*>();
  }
  fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  fastjet::contrib::SoftDrop sd(m_SDBeta, m_SDZCut);
  if (m_Verbosity > 5)
  {
    std::cout << "FastJetAlgo::get_jets : created SoftDrop groomer configuration : " << sd.description() << std::endl;
  }

  // translate into jet output...
  std::vector<Jet*> jets;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    Jet* jet = new Jetv1();
    jet->set_px(fastjets[ijet].px());
    jet->set_py(fastjets[ijet].py());
    jet->set_pz(fastjets[ijet].pz());
    jet->set_e(fastjets[ijet].e());
    jet->set_id(ijet);

    // if SoftDrop enabled, and jets have > 5 GeV (do not waste time
    // on very low-pT jets), run SD and pack output into jet properties
    if (m_SDFlag && fastjets[ijet].perp() > 5)
    {
      fastjet::PseudoJet sd_jet = sd(fastjets[ijet]);

      if (m_Verbosity > 5)
      {
        std::cout << "original    jet: pt / eta / phi / m = " << fastjets[ijet].perp() << " / " << fastjets[ijet].eta() << " / " << fastjets[ijet].phi() << " / " << fastjets[ijet].m() << std::endl;
        std::cout << "SoftDropped jet: pt / eta / phi / m = " << sd_jet.perp() << " / " << sd_jet.eta() << " / " << sd_jet.phi() << " / " << sd_jet.m() << std::endl;

        std::cout << "  delta_R between subjets: " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() << std::endl;
        std::cout << "  symmetry measure(z):     " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry() << std::endl;
        std::cout << "  mass drop(mu):           " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu() << std::endl;
      }

      // attach SoftDrop quantities as jet properties
      jet->set_property(Jet::PROPERTY::prop_zg, sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
      jet->set_property(Jet::PROPERTY::prop_Rg, sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R());
      jet->set_property(Jet::PROPERTY::prop_mu, sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu());
    }

    // copy components into output jet
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (auto & comp : comps)
    {
      Jet* particle = particles[comp.user_index()];

      for (Jet::Iter iter = particle->begin_comp();
           iter != particle->end_comp();
           ++iter)
      {
        jet->insert_comp(iter->first, iter->second);
      }
    }

    jets.push_back(jet);
  }

  if (m_Verbosity > 1) std::cout << "FastJetAlgo::process_event -- exited" << std::endl;

  return jets;
}

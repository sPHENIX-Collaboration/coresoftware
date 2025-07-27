
#include "FastJetAlgoSub.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

// standard includes
#include <iostream>
#include <sstream>  // for basic_ostringstream
#include <vector>

FastJetAlgoSub::FastJetAlgoSub(const FastJetOptions& options)
  : m_opt{options}
{
  fastjet::ClusterSequence clusseq;
  if (m_opt.verbosity > 0)
  {
    fastjet::ClusterSequence::print_banner();
  }
  else
  {
    std::ostringstream nullstream;
    fastjet::ClusterSequence::set_fastjet_banner_stream(&nullstream);
    fastjet::ClusterSequence::print_banner();
    fastjet::ClusterSequence::set_fastjet_banner_stream(&std::cout);
  }
}

void FastJetAlgoSub::identify(std::ostream& os)
{
  os << "   FastJetAlgoSub: ";
  if (m_opt.algo == Jet::ANTIKT)
  {
    os << "ANTIKT r=" << m_opt.jet_R;
  }
  else if (m_opt.algo == Jet::KT)
  {
    os << "KT r=" << m_opt.jet_R;
  }
  else if (m_opt.algo == Jet::CAMBRIDGE)
  {
    os << "CAMBRIDGE r=" << m_opt.jet_R;
  }
  os << std::endl;
}

/* std::vector<Jet*> FastJetAlgoSub::get_jets(std::vector<Jet*> particles) */
/* { }; //  deprecated by iterating from JetMap; now is JetContainer and most code moved into */
//  into cluster_and_fill

void FastJetAlgoSub::cluster_and_fill(std::vector<Jet*>& particles, JetContainer* jetcont)
{
  if (m_opt.verbosity > 1)
  {
    std::cout << "FastJetAlgoSub::process_event -- entered" << std::endl;
  }

  // translate to fastjet
  std::vector<fastjet::PseudoJet> pseudojets;
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart)
  {
    float this_e = particles[ipart]->get_e();

    if (this_e == 0.)
    {
      continue;
    }

    float this_px = particles[ipart]->get_px();
    float this_py = particles[ipart]->get_py();
    float this_pz = particles[ipart]->get_pz();

    if (this_e < 0)
    {
      // make energy = +1 MeV for purposes of clustering
      float e_ratio = 0.001 / this_e;

      this_e = this_e * e_ratio;
      this_px = this_px * e_ratio;
      this_py = this_py * e_ratio;
      this_pz = this_pz * e_ratio;

      if (m_opt.verbosity > 5)
      {
        std::cout << " FastJetAlgoSub input particle with negative-E, original kinematics px / py / pz / E = ";
        std::cout << particles[ipart]->get_px() << " / " << particles[ipart]->get_py() << " / " << particles[ipart]->get_pz() << " / " << particles[ipart]->get_e() << std::endl;
        std::cout << " --> entering with modified kinematics px / py / pz / E = " << this_px << " / " << this_py << " / " << this_pz << " / " << this_e << std::endl;
      }
    }

    fastjet::PseudoJet pseudojet(this_px, this_py, this_pz, this_e);

    pseudojet.set_user_index(ipart);
    pseudojets.push_back(pseudojet);
  }

  // run fast jet
  fastjet::JetDefinition* jetdef = nullptr;
  if (m_opt.algo == Jet::ANTIKT)
  {
    jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_opt.algo == Jet::KT)
  {
    jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_opt.algo == Jet::CAMBRIDGE)
  {
    jetdef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else
  {
    return;
  }

  fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  // translate into jet output...
  std::vector<Jet*> jets;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    auto* jet = jetcont->add_jet();

    if (m_opt.verbosity > 5 && fastjets[ijet].perp() > 15)
    {
      std::cout << " FastJetAlgoSub : jet # " << ijet << " comes out of clustering with pt / eta / phi = " << fastjets[ijet].perp() << " / " << fastjets[ijet].eta() << " / " << fastjets[ijet].phi();
      std::cout << ", px / py / pz / e = " << fastjets[ijet].px() << " / " << fastjets[ijet].py() << " / " << fastjets[ijet].pz() << " / " << fastjets[ijet].e() << std::endl;
    }

    float total_px = 0;
    float total_py = 0;
    float total_pz = 0;
    float total_e = 0;

    // copy components into output jet
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (auto& comp : comps)
    {
      Jet* particle = particles[comp.user_index()];

      total_px += particle->get_px();
      total_py += particle->get_py();
      total_pz += particle->get_pz();
      total_e += particle->get_e();
      jet->insert_comp(particle->get_comp_vec(), true);
    }

    jet->set_comp_sort_flag();  // make sure jet know comps might not be sorted
                                // alternatively can just ommit the `true`
                                // in insert_comp call above
    jet->set_px(total_px);
    jet->set_py(total_py);
    jet->set_pz(total_pz);
    jet->set_e(total_e);
    jet->set_id(ijet);

    if (m_opt.verbosity > 5 && fastjets[ijet].perp() > 15)
    {
      std::cout << " FastJetAlgoSub : jet # " << ijet << " after correcting for proper constituent kinematics, pt / eta / phi = " << jet->get_pt() << " / " << jet->get_eta() << " / " << jet->get_phi();
      std::cout << ", px / py / pz / e = " << jet->get_px() << " / " << jet->get_py() << " / " << jet->get_pz() << " / " << jet->get_e() << std::endl;
    }
  }

  if (m_opt.verbosity > 1)
  {
    std::cout << "FastJetAlgoSub::process_event -- exited" << std::endl;
  }
  return;
}

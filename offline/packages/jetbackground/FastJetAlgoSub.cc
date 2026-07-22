
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

namespace
{
  // which calorimeter layer a jet input source belongs to,
  // for the EMCal/iHCal/oHCal jet energy fractions
  enum class CaloLayer
  {
    NONE,
    EMCAL,
    IHCAL,
    OHCAL
  };

  CaloLayer get_calo_layer(Jet::SRC src)
  {
    switch (src)
    {
    case Jet::CEMC_TOWER:
    case Jet::CEMC_CLUSTER:
    case Jet::CEMC_TOWER_RETOWER:
    case Jet::CEMC_TOWER_SUB1:
    case Jet::CEMC_TOWER_SUB1CS:
    case Jet::CEMC_TOWERINFO:
    case Jet::CEMC_TOWERINFO_RETOWER:
    case Jet::CEMC_TOWERINFO_SUB1:
    case Jet::CEMC_TOWERINFO_EMBED:
    case Jet::CEMC_TOWERINFO_SIM:
    case Jet::ECAL_TOPO_CLUSTER:
      return CaloLayer::EMCAL;
    case Jet::HCALIN_TOWER:
    case Jet::HCALIN_CLUSTER:
    case Jet::HCALIN_TOWER_SUB1:
    case Jet::HCALIN_TOWER_SUB1CS:
    case Jet::HCALIN_TOWERINFO:
    case Jet::HCALIN_TOWERINFO_SUB1:
    case Jet::HCALIN_TOWERINFO_EMBED:
    case Jet::HCALIN_TOWERINFO_SIM:
      return CaloLayer::IHCAL;
    case Jet::HCALOUT_TOWER:
    case Jet::HCALOUT_CLUSTER:
    case Jet::HCALOUT_TOWER_SUB1:
    case Jet::HCALOUT_TOWER_SUB1CS:
    case Jet::HCALOUT_TOWERINFO:
    case Jet::HCALOUT_TOWERINFO_SUB1:
    case Jet::HCALOUT_TOWERINFO_EMBED:
    case Jet::HCALOUT_TOWERINFO_SIM:
      return CaloLayer::OHCAL;
    default:
      return CaloLayer::NONE;
    }
  }
}  // namespace

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
    float w_t_sum = 0;
    float w_e_sum = 0;
    float emcal_e = 0;
    float ihcal_e = 0;
    float ohcal_e = 0;
    // copy components into output jet
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (auto& comp : comps)
    {
      Jet* particle = particles[comp.user_index()];

      total_px += particle->get_px();
      total_py += particle->get_py();
      total_pz += particle->get_pz();
      total_e += particle->get_e();
      if (m_opt.calc_calo_fracs && !particle->get_comp_vec().empty())
      {
        switch (get_calo_layer(particle->get_comp_vec().front().first))
        {
        case CaloLayer::EMCAL:
          emcal_e += particle->get_e();
          break;
        case CaloLayer::IHCAL:
          ihcal_e += particle->get_e();
          break;
        case CaloLayer::OHCAL:
          ohcal_e += particle->get_e();
          break;
        case CaloLayer::NONE:
          break;
        }
      }
      if(particle->size_properties() > Jet::PROPERTY::prop_t)
	{
	  if(!std::isnan(particle->get_property(Jet::PROPERTY::prop_t)))
	    {
	      w_t_sum += particle->get_property(Jet::PROPERTY::prop_t) * particle->get_e();
	      w_e_sum += particle->get_e();
	    }
	}
	jet->insert_comp(particle->get_comp_vec(), true);
    }
    if(jet->size_properties() < Jet::PROPERTY::prop_t+1)
      {
	jet->resize_properties(Jet::PROPERTY::prop_t + 1);
      }
    jet->set_property(Jet::PROPERTY::prop_t,w_t_sum/w_e_sum); //This intentionally becomes nan (0/0) if the
                                                              //value has not been filled at all so that
                                                              //the jet auto-fails timing cuts if necessary

    jet->set_comp_sort_flag();  // make sure jet know comps might not be sorted
                                // alternatively can just ommit the `true`
                                // in insert_comp call above
    jet->set_px(total_px);
    jet->set_py(total_py);
    jet->set_pz(total_pz);
    jet->set_e(total_e);
    jet->set_id(ijet);

    // calo energy fractions: original tower energy per layer / jet total energy
    if (m_opt.calc_calo_fracs && jet->get_e() != 0)
    {
      jet->set_emcal_frac(emcal_e / jet->get_e());
      jet->set_ihcal_frac(ihcal_e / jet->get_e());
      jet->set_ohcal_frac(ohcal_e / jet->get_e());
    }

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

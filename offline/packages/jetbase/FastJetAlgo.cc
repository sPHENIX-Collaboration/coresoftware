#include "FastJetAlgo.h"

#include "Jet.h"
#include "Jetv2.h"
#include "JetContainer.h"

#include <phool/phool.h>

#include <TSystem.h>
#include <TClonesArray.h>

// fastjet includes
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>  // for FunctionOfPse...
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/RecursiveSymmetryCutBase.hh>  // for RecursiveSymm...
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/contrib/ConstituentSubtractor.hh>


// standard includes
#include <cmath>  // for isfinite
#include <iostream>
#include <map>      // for _Rb_tree_iterator
#include <memory>   // for allocator_traits<>::value_type
#include <string>   // for operator<<
#include <utility>  // for pair
#include <vector>
#include <fstream>
#include <cassert>

FastJetAlgo::FastJetAlgo(const FastJetOptions& options) :
  m_opt { options }
{
  fastjet::ClusterSequence clusseq;
  if (m_opt.verbosity > 0)
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

  m_opt.print(os);
}

fastjet::JetDefinition FastJetAlgo::get_fastjet_definition() {
  if (m_opt.algo == Jet::ANTIKT)
  {
    return fastjet::JetDefinition(fastjet::antikt_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_opt.algo == Jet::KT)
  {
    return fastjet::JetDefinition(fastjet::kt_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_opt.algo == Jet::CAMBRIDGE)
  {
    return fastjet::JetDefinition(fastjet::cambridge_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
  else
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "Warning, no recognized jet clustering algorithm provided in FastJetAlgo" << std::endl
              << "defaulting to antikt_algorithm" << std::endl;
    //return a dummy definition
    return fastjet::JetDefinition(fastjet::antikt_algorithm, m_opt.jet_R, fastjet::E_scheme, fastjet::Best);
  }
}

fastjet::Selector FastJetAlgo::get_selector() {
  // only selectors available are jet_min_pt and jet_max_eta
  if (m_opt.use_jet_max_eta && m_opt.use_jet_min_pt) {
    return fastjet::SelectorAbsRapMax(m_opt.jet_max_eta) 
       &&  fastjet::SelectorPtMin(m_opt.jet_min_pt);
  } else if (m_opt.use_jet_max_eta) {
    return fastjet::SelectorAbsRapMax(m_opt.jet_max_eta);
  } else {
    return fastjet::SelectorPtMin(m_opt.jet_min_pt);
  }
}

std::vector<fastjet::PseudoJet> FastJetAlgo::cluster_jets(
    std::vector<fastjet::PseudoJet>& pseudojets
) {
  auto jetdef = get_fastjet_definition();
  m_cluseq = new fastjet::ClusterSequence( pseudojets, jetdef );

  if (m_opt.use_jet_selection) {
    auto selector = get_selector();
    return fastjet::sorted_by_pt(selector(m_cluseq->inclusive_jets()));
  } else {
    return fastjet::sorted_by_pt(m_cluseq->inclusive_jets());
  }
}

std::vector<fastjet::PseudoJet> FastJetAlgo::cluster_area_jets(
    std::vector<fastjet::PseudoJet>& pseudojets
) {

  auto jetdef = get_fastjet_definition();

  fastjet::AreaDefinition area_def ( 
        fastjet::active_area_explicit_ghosts, 
        fastjet::GhostedAreaSpec(m_opt.ghost_max_rap, 1, m_opt.ghost_area)
  );

  m_cluseqarea = new fastjet::ClusterSequenceArea( pseudojets, jetdef, area_def );

  fastjet::Selector selector = (
      m_opt.use_jet_selection 
      ? (!fastjet::SelectorIsPureGhost() && get_selector())
      : !fastjet::SelectorIsPureGhost()
  );

  return fastjet::sorted_by_pt(selector(m_cluseqarea->inclusive_jets()));

}

float FastJetAlgo::calc_rhomeddens(std::vector<fastjet::PseudoJet>& constituents) {
  fastjet::AreaDefinition area_def ( 
      fastjet::active_area_explicit_ghosts, 
      fastjet::GhostedAreaSpec(m_opt.ghost_max_rap, 1, m_opt.ghost_area)
      );

  fastjet::Selector rho_select =  (!fastjet::SelectorNHardest(m_opt.nhardestcut_jetmedbkgdens)) 
    * fastjet::SelectorAbsEtaMax(m_opt.etahardestcut_jetmedbkgdens); // <--

  fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, m_opt.jet_R); // <--
  fastjet::JetMedianBackgroundEstimator bge {rho_select, jet_def_bkgd, area_def};
  bge.set_particles(constituents);
  return bge.rho();
}

std::vector<fastjet::PseudoJet> 
    FastJetAlgo::jets_to_pseudojets(std::vector<Jet*>& particles) {
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
    if (m_opt.use_constituent_min_pt && pseudojet.perp() < m_opt.constituent_min_pt) continue;
    pseudojet.set_user_index(ipart);
    pseudojets.push_back(pseudojet);
  }
  return pseudojets;
}

void FastJetAlgo::first_call_init(JetContainer* jetcont) {
  m_first_cluster_call = false;
  m_opt.initialize();

  if (jetcont == nullptr) return;

  if (m_opt.doSoftDrop) {
    jetcont->add_property( {Jet::PROPERTY::prop_zg, Jet::PROPERTY::prop_Rg, Jet::PROPERTY::prop_mu} );
    m_zg_index = jetcont->property_index(Jet::PROPERTY::prop_zg);
    m_Rg_index = jetcont->property_index(Jet::PROPERTY::prop_Rg);
    m_mu_index = jetcont->property_index(Jet::PROPERTY::prop_mu);
  }

  if (m_opt.calc_area) {
    jetcont->add_property(Jet::PROPERTY::prop_area);
    m_area_index = jetcont->property_index(Jet::PROPERTY::prop_area);
  }

  if (m_opt.cs_calc_constsub) {
    cs_bge_rho =  new fastjet::GridMedianBackgroundEstimator(m_opt.cs_max_eta, 0.5);
    cs_subtractor = new fastjet::contrib::ConstituentSubtractor();
	  cs_subtractor->set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); 
	  cs_subtractor->set_max_distance(m_opt.cs_max_dist); // free parameter for the maximal allowed distance between particle i and ghost k
	  cs_subtractor->set_alpha(m_opt.cs_alpha);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
	  cs_subtractor->set_ghost_area(m_opt.cs_ghost_area); // free parameter for the density of ghosts. 
	  cs_subtractor->set_max_eta(m_opt.cs_max_eta); // parameter for the maximal eta cut
	  cs_subtractor->set_background_estimator(cs_bge_rho); // specify the background estimator to estimate rho.
    if (m_opt.cs_max_pt > 0) {
      cs_sel_max_pt = new fastjet::Selector(fastjet::SelectorPtMax(m_opt.cs_max_pt));
      cs_subtractor->set_particle_selector(cs_sel_max_pt);
    }
    cs_subtractor->initialize();
    if (m_opt.verbosity > 1) std::cout << cs_subtractor->description() << std::endl;
  }

  jetcont->set_algo(m_opt.algo);
  jetcont->set_jetpar_R(m_opt.jet_R);
}

void FastJetAlgo::cluster_and_fill(std::vector<Jet*>& particles, JetContainer* jetcont)
{
  if (m_first_cluster_call) first_call_init(jetcont);
    // initalize the properties in JetContainer

  if (m_opt.verbosity > 1) std::cout << "   Verbosity>1 FastJetAlgo::process_event -- entered" << std::endl;
  if (m_opt.verbosity > 8) std::cout << "   Verbosity>8 #input particles: " << particles.size() << std::endl;

  // translate input jets to input fastjets
  auto pseudojets = jets_to_pseudojets(particles);

  // if using constituent subtraction, oberve maximum eta and subtract the constituents
  if (m_opt.cs_calc_constsub) {
    if (m_opt.verbosity > 100) {
      std::cout << " Before Constituent Subtraction: " << std::endl;
      int i = 0;
      double sumpt = 0.;
      for (auto c : pseudojets) {
        sumpt += c.perp();
        if (i<100) std::cout << Form(" jet[%2i] %8.4f  sum %8.4f", i, c.perp(), sumpt) << std::endl;
        i++;
      }
      auto _c = pseudojets.back();
      std::cout << Form(" jet[%2i] %8.4f  sum %8.4f", i++, _c.perp(), sumpt) << std::endl << std::endl;
    }

    pseudojets = fastjet::SelectorAbsEtaMax(m_opt.cs_max_eta)(pseudojets);
    cs_bge_rho->set_particles(pseudojets);
    auto subtracted_pseudojets = cs_subtractor->subtract_event(pseudojets);
    pseudojets = std::move(subtracted_pseudojets);
    
    if (m_opt.verbosity > 100) {
      std::cout << " After Constituent Subtraction: " << std::endl;
      int i = 0;
      double sumpt = 0.;
      for (auto c : pseudojets) {
        sumpt += c.perp();
        if (i<100) std::cout << Form(" jet[%2i] %8.4f  sum %8.4f", i, c.perp(), sumpt) << std::endl;
        i++;
      }
      auto _c = pseudojets.back();
      std::cout << Form(" jet[%2i] %8.4f  sum %8.4f", i++, _c.perp(), sumpt) << std::endl << std::endl;
    }
  }

  if (m_opt.calc_jetmedbkgdens) jetcont->set_rho_median(calc_rhomeddens(pseudojets));

  auto fastjets = (m_opt.calc_area ? cluster_area_jets(pseudojets) : cluster_jets(pseudojets));

  if (m_opt.verbosity > 8) std::cout << "   Verbosity>8 fastjets: " << fastjets.size() << std::endl;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    auto* jet = jetcont->add_jet(); // put a new Jetv2 into the TClonesArray
    jet->set_px(fastjets[ijet].px());
    jet->set_py(fastjets[ijet].py());
    jet->set_pz(fastjets[ijet].pz());
    jet->set_e(fastjets[ijet].e());
    jet->set_id(ijet);

    if (m_opt.calc_area) {
      jet->set_property(m_area_index, fastjets[ijet].area());
    }

    // if SoftDrop enabled, and jets have > 5 GeV (do not waste time
    // on very low-pT jets), run SD and pack output into jet properties
    if (m_opt.doSoftDrop && fastjets[ijet].perp() > 5)
    {
      fastjet::contrib::SoftDrop sd(m_opt.SD_beta, m_opt.SD_zcut);
      if (m_opt.verbosity > 5)
      {
        std::cout << "FastJetAlgo::get_jets : created SoftDrop groomer configuration : " 
          << sd.description() << std::endl;
      }

      fastjet::PseudoJet sd_jet = sd(fastjets[ijet]);

      if (m_opt.verbosity > 5)
      {
        std::cout << "original    jet: pt / eta / phi / m = " << fastjets[ijet].perp() 
            << " / " << fastjets[ijet].eta() << " / " << fastjets[ijet].phi() << " / " 
            << fastjets[ijet].m() << std::endl;
        std::cout << "SoftDropped jet: pt / eta / phi / m = " << sd_jet.perp() << " / " 
            << sd_jet.eta() << " / " << sd_jet.phi() << " / " << sd_jet.m() << std::endl;

        std::cout << "  delta_R between subjets: " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() << std::endl;
        std::cout << "  symmetry measure(z):     " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry() << std::endl;
        std::cout << "  mass drop(mu):           " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu() << std::endl;
      }
 
      // attach SoftDrop quantities as jet properties
      jet->set_property(m_zg_index, sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
      jet->set_property(m_Rg_index, sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() );
      jet->set_property(m_mu_index, sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu()      );
    }
 
    // Count clustered components. If desired, put original components into the output jet.
    int n_clustered = 0;
    std::vector<fastjet::PseudoJet> constituents = fastjets[ijet].constituents();
    if (m_opt.calc_area) {
      for (auto& comp : constituents) {
        if (comp.is_pure_ghost()) continue;
        ++n_clustered;
        if (m_opt.save_jet_components) {
          jet->insert_comp(particles[comp.user_index()]->get_comp_vec(), true);
        }
      } // end loop over all constituents
    } else { // didn't calculate jet area
      n_clustered += constituents.size(); 
      if (m_opt.save_jet_components) {
        for (auto& comp : constituents) {
          jet->insert_comp(particles[comp.user_index()]->get_comp_vec(), true);
        }
      }
    }
    jet->set_comp_sort_flag(); // make surce comp knows it might not be sorted
  }
  if (m_opt.verbosity > 1) std::cout << "FastJetAlgo::process_event -- exited" << std::endl;
  delete (m_opt.calc_area ? m_cluseqarea : m_cluseq); //if (m_cluseq) delete m_cluseq;
}


std::vector<Jet*> FastJetAlgo::get_jets(std::vector<Jet*> particles)
{
  // translate to fastjet
  auto pseudojets = jets_to_pseudojets(particles);
  auto fastjets = cluster_jets(pseudojets);

  fastjet::contrib::SoftDrop sd(m_opt.SD_beta, m_opt.SD_zcut);
  if (m_opt.verbosity > 5)
  {
    std::cout << "FastJetAlgo::get_jets : created SoftDrop groomer configuration : " << sd.description() << std::endl;
  }

  // translate into jet output...
  std::vector<Jet*> jets;
  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    Jet* jet = new Jetv2();
    jet->set_px(fastjets[ijet].px());
    jet->set_py(fastjets[ijet].py());
    jet->set_pz(fastjets[ijet].pz());
    jet->set_e(fastjets[ijet].e());
    jet->set_id(ijet);

    // if SoftDrop enabled, and jets have > 5 GeV (do not waste time
    // on very low-pT jets), run SD and pack output into jet properties
    if (m_opt.doSoftDrop && fastjets[ijet].perp() > 5)
    {
      fastjet::PseudoJet sd_jet = sd(fastjets[ijet]);

      if (m_opt.verbosity > 5)
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

    // Count clustered components. If desired, put original components into the output jet.
    int n_clustered = 0;
    std::vector<fastjet::PseudoJet> constituents = fastjets[ijet].constituents();
    if (m_opt.calc_area) {
      for (auto& comp : constituents) {
        if (comp.is_pure_ghost()) continue;
        ++n_clustered;
        if (m_opt.save_jet_components) {
          jet->insert_comp(particles[comp.user_index()]->get_comp_vec(), true);
        }
      } // end loop over all constituents
    } else { // didn't save jet area
      n_clustered += constituents.size(); 
      if (m_opt.save_jet_components) {
        for (auto& comp : constituents) {
          jet->insert_comp(particles[comp.user_index()]->get_comp_vec(), true);
        }
      }
    }
    /* jet->set_n_clustered(n_clustered); */
    jet->set_comp_sort_flag(); // make surce comp knows it might not be sorted
                               // can alternatively just remove the `true` parameter
                               // from the insert_comp function calls above
    jets.push_back(jet);
  }

  /* std::sort(jets.begin(), jets.end(), [](Jet*a, Jet*b){ return a->get_et() > b->get_et();} ); */

  if (m_opt.verbosity > 1) std::cout << "FastJetAlgo::process_event -- exited" << std::endl;

  delete m_cluseq;

  return jets;
}


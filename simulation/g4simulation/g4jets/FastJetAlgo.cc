
#include "FastJetAlgo.h"

#include "Jet.h"
#include "Jetv1.h"

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

// SoftDrop includes
#include <fastjet/contrib/SoftDrop.hh>

// standard includes
#include <iostream>
#include <map>                         // for _Rb_tree_iterator
#include <memory>                      // for allocator_traits<>::value_type
#include <utility>                     // for pair
#include <vector>

FastJetAlgo::FastJetAlgo(Jet::ALGO algo, float par, int verbosity)
  : _verbosity(verbosity)
  , _algo(algo)
  , _par(par)
  , _do_SD(false)
  , _SD_beta(0)
  , _SD_zcut(0.1)
{
  fastjet::ClusterSequence clusseq;
  if (_verbosity > 0)
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
  if (_algo == Jet::ANTIKT)
    os << "ANTIKT r=" << _par;
  else if (_algo == Jet::KT)
    os << "KT r=" << _par;
  else if (_algo == Jet::CAMBRIDGE)
    os << "CAMBRIDGE r=" << _par;
  os << std::endl;
}

std::vector<Jet*> FastJetAlgo::get_jets(std::vector<Jet*> particles)
{
  if (_verbosity > 1) std::cout << "FastJetAlgo::process_event -- entered" << std::endl;

  // translate to fastjet
  std::vector<fastjet::PseudoJet> pseudojets;
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart)
  {
    // fastjet performs strangely with exactly (px,py,pz,E) =
    // (0,0,0,0) inputs, such as placeholder towers or those with
    // zero'd out energy after CS. this catch also in FastJetAlgoSub
    if (particles[ipart]->get_e() == 0.) continue;

    fastjet::PseudoJet pseudojet(particles[ipart]->get_px(),
                                 particles[ipart]->get_py(),
                                 particles[ipart]->get_pz(),
                                 particles[ipart]->get_e());
    pseudojet.set_user_index(ipart);
    pseudojets.push_back(pseudojet);
  }

  // run fast jet
  fastjet::JetDefinition* jetdef = nullptr;
  if (_algo == Jet::ANTIKT)
    jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm, _par, fastjet::E_scheme, fastjet::Best);
  else if (_algo == Jet::KT)
    jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, _par, fastjet::E_scheme, fastjet::Best);
  else if (_algo == Jet::CAMBRIDGE)
    jetdef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, _par, fastjet::E_scheme, fastjet::Best);
  else
    return std::vector<Jet*>();
  fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();
  delete jetdef;

  fastjet::contrib::SoftDrop sd( _SD_beta, _SD_zcut );
  if ( _verbosity > 5 )
    std::cout << "FastJetAlgo::get_jets : created SoftDrop groomer configuration : " << sd.description() << std::endl;

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
    if ( _do_SD && fastjets[ijet].perp() > 5 ) {

      fastjet::PseudoJet sd_jet = sd( fastjets[ijet] );
            
      if ( _verbosity > 5 ) {
	std::cout << "original    jet: pt / eta / phi / m = " << fastjets[ijet].perp() << " / " <<  fastjets[ijet].eta() << " / " << fastjets[ijet].phi() << " / " << fastjets[ijet].m() << std::endl;
	std::cout << "SoftDropped jet: pt / eta / phi / m = " << sd_jet.perp() << " / " <<  sd_jet.eta() << " / " << sd_jet.phi() << " / " << sd_jet.m() << std::endl;
	
	std::cout << "  delta_R between subjets: " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() << std::endl;
	std::cout << "  symmetry measure(z):     " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry() << std::endl;
	std::cout << "  mass drop(mu):           " << sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu() << std::endl;      
      }

      // attach SoftDrop quantities as jet properties
      jet->set_property( Jet::PROPERTY::prop_zg , sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry() );
      jet->set_property( Jet::PROPERTY::prop_Rg , sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() );
      jet->set_property( Jet::PROPERTY::prop_mu , sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu() );

    }

    // copy components into output jet
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (unsigned int icomp = 0; icomp < comps.size(); ++icomp)
    {
      Jet* particle = particles[comps[icomp].user_index()];

      for (Jet::Iter iter = particle->begin_comp();
           iter != particle->end_comp();
           ++iter)
      {
        jet->insert_comp(iter->first, iter->second);
      }
    }

    jets.push_back(jet);
  }

  if (_verbosity > 1) std::cout << "FastJetAlgo::process_event -- exited" << std::endl;

  return jets;
}

#ifndef JETBASE_FASTJETALGO_H
#define JETBASE_FASTJETALGO_H

#include "Jet.h"
#include "JetAlgo.h"
#include "FastJetOptions.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

#include <cmath>     // for NAN
#include <climits>     // for NAN
#include <iostream>  // for cout, ostream
#include <vector>    // for vector
                     

namespace fastjet {
  class PseudoJet;
  class GridMedianBackgroundEstimator;
  class SelectorPtMax;
  namespace contrib {
    class ConstituentSubtractor;
  }
}

class JetContainer;

class FastJetAlgo : public JetAlgo
{
 public:
  FastJetAlgo(const FastJetOptions& options);
  ~FastJetAlgo() override {}

  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return m_opt.algo; }

  //----------------------------------------------------------------------
  //  Legacy code interface. It is better to use FastJetOptions, but
  //  there is no harm is using these, as well.
  //----------------------------------------------------------------------
  FastJetAlgo(Jet::ALGO algo, float par, int verbosity=0) :
     FastJetAlgo({{algo, JET_R,par,VERBOSITY, static_cast<float>(verbosity)}})
  {}
  void set_do_SoftDrop   (bool  do_SD) { m_opt.doSoftDrop = do_SD; }
  void set_SoftDrop_beta (float beta)  { m_opt.SD_beta    = beta;  }
  void set_SoftDrop_zcut (float zcut)  { m_opt.SD_zcut    = zcut;  }
  //--end-legacy-code-interface-------------------------------------------

  std::vector<Jet*> get_jets(std::vector<Jet*> particles) override;
  void cluster_and_fill(std::vector<Jet*>& part_in, JetContainer* jets_out) override;

 private:
  FastJetOptions m_opt {};
  bool m_first_cluster_call { true };

  // for convenience save indices of the zg, Rg, mu, and area for jets in the JetContainer
  Jet::PROPERTY m_zg_index   { Jet::PROPERTY::no_property };
  Jet::PROPERTY m_Rg_index   { Jet::PROPERTY::no_property };
  Jet::PROPERTY m_mu_index   { Jet::PROPERTY::no_property };
  Jet::PROPERTY m_area_index { Jet::PROPERTY::no_property };

  // Internal processes
  std::vector<fastjet::PseudoJet> jets_to_pseudojets(std::vector<Jet*>& particles);
  std::vector<fastjet::PseudoJet> cluster_jets(std::vector<fastjet::PseudoJet>& constituents);
  std::vector<fastjet::PseudoJet> cluster_area_jets(std::vector<fastjet::PseudoJet>& constituents);
  float calc_rhomeddens(std::vector<fastjet::PseudoJet>& constituents);
  fastjet::JetDefinition get_fastjet_definition();
  fastjet::Selector get_selector();
  void first_call_init(JetContainer*_=nullptr);

  // private members
  fastjet::contrib::ConstituentSubtractor* cs_subtractor = nullptr;
  fastjet::GridMedianBackgroundEstimator*  cs_bge_rho    = nullptr;
  fastjet::Selector*                       cs_sel_max_pt = nullptr;


  fastjet::ClusterSequence* m_cluseq {nullptr};
  fastjet::ClusterSequence* m_cluseqarea {nullptr};
};

#endif

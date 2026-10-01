/*!
 * \file DetermineTowerRho.h
 * \brief UE background rho calculator.
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 2.0.1 $
 * \date $Date: 02/01/2024. Revised 09/19/2024$
 */

#ifndef JETBACKGROUND_DETERMINETOWERHO_H
#define JETBACKGROUND_DETERMINETOWERHO_H

#include "TowerRho.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class JetContainer;
class JetInput;
class Jet;

namespace fastjet
{
  class PseudoJet;
  class Selector;
}  // namespace fastjet

class DetermineTowerRho : public SubsysReco
{
 public:
  DetermineTowerRho(const std::string &name = "DetermineTowerRho");
  ~DetermineTowerRho() override;

  // standard Fun4All methods
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  // add rho method (Area or Multiplicity)
  //
  // The background jets are clustered with fastjet and selected with the jet
  // acceptance (set_jet_abs_eta). The surviving jets are then converted back
  // to sPHENIX Jets: the four-vector of each is rebuilt from the ORIGINAL
  // constituents,
  // Every selected jet -- the omitted hardest ones (the seeds) included -- is
  // written to the node tree as a JetContainer named jet_node. The default rho
  // node is TowerRho_<METHOD> and the default jet node is <rho node>_SeedJets,
  // e.g. TowerRho_MULT -> TowerRho_MULT_SeedJets, or with output_node =
  // "TowerRho_CEMC_MULT" -> TowerRho_CEMC_MULT_SeedJets.
  //   Jet::PROPERTY::prop_SeedItr : 1 if it was omitted as one of the hardest
  //                                 (a seed), 0 if it entered rho
  //   Jet::PROPERTY::prop_signedSumeT : the signed scalar sum of its constituent
  //                                 pT, the quantity rho is estimated from

  void add_method(
    TowerRho::Method rho_method,
    const std::string &output_node = "",
    const std::string &jet_node = ""
  );

  // inputs for estimating background
  void add_input(JetInput *input) { m_inputs.push_back(input); }
  void add_tower_input(JetInput *input) { add_input(input); }  // for backwards compatibility

  // -- background jet clustering --
  // algorithm (default KT) and radius parameter R (default 0.4)
  void set_algo(const Jet::ALGO algo) { m_bkgd_jet_algo = algo; }
  Jet::ALGO get_algo() const { return m_bkgd_jet_algo; }

  void set_par(const float val) { m_par = val; }
  float get_par() const { return m_par; }

  // -- ghosts (AREA method only) --
  // ghost area (default 0.01)
  void set_ghost_area(const float val) { m_ghost_area = val; }
  float get_ghost_area() const { return m_ghost_area; }

  // -- background jet acceptance --
  void set_jet_abs_eta(const float abseta) { m_abs_jet_eta_range = abseta; }
  float get_jet_abs_eta() const { return m_abs_jet_eta_range; }

  // number of hardest jets, ranked by prop_signedSumeT, left out of the rho
  // estimate (default 2)
  void set_omit_nhardest(const unsigned int val) { m_omit_nhardest = val; }
  unsigned int get_omit_nhardest() const { return m_omit_nhardest; }

 private:

  // variables
  std::vector<JetInput *> m_inputs{};

  std::vector<std::string> m_output_nodes{};
  std::vector<std::string> m_jet_output_nodes{};

  std::vector<TowerRho::Method> m_rho_methods{};

  Jet::ALGO m_bkgd_jet_algo{Jet::ALGO::KT};
  float m_par{0.4};

  float m_ghost_area{0.01};

  float m_abs_jet_eta_range{-1.0};  // negative = event tower eta range shrunk by R
  unsigned int m_omit_nhardest{2};

  // internal methods
  int CreateNodes(PHCompositeNode *topNode);

  // Convert the selected fastjet jets back to sPHENIX Jets
  void ConvertJets(
    JetContainer *jets,
    const std::vector<fastjet::PseudoJet> &fastjets,
    const std::vector<Jet *> &particles,
    bool with_area
  ) const;


  // prop_signedSumeT of jet ijet: the per-jet quantity the median is taken over,
  static float SignedSumeT(JetContainer *jets, unsigned int ijet);

  // Drop the m_omit_nhardest hardest sPHENIX jets
  std::vector<unsigned int> SelectJets(JetContainer *jets) const;

  // Estimate rho and sigma from the sPHENIX jets of `jets` listed in `keep`
  void CalcRho(
    JetContainer *jets,
    const std::vector<unsigned int> &keep,
    TowerRho::Method rho_method,
    float &rho,
    float &sigma
  ) const;

  static float CalcPercentile(
    const std::vector<float> &sorted_vec,
    const float percentile,
    const float nempty
  );

  static void CalcMedianStd(
    const std::vector<float> &vec,
    float n_empty_jets,
    float &median,
    float &std_dev
  );

  fastjet::Selector get_jet_selector(double eta_min, double eta_max) const;
};

#endif  // JETBACKGROUND_DETERMINETOWERHO_H

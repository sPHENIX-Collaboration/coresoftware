#ifndef RHOBASE_DETERMINETOWERRHO_H
#define RHOBASE_DETERMINETOWERRHO_H

//===========================================================
/// \file DetermineTowerRho.h
/// \brief UE background rho calculator
/// \author Tanner Mengel
//===========================================================

#include "TowerRho.h"
#include "TowerRhov1.h"

#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetInput.h>

// fun4all includes
#include <fun4all/SubsysReco.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

// system includes
#include <ostream>
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;

/// \class DetermineTowerRho
///
/// \brief UE background calculator
///
/// This module estimates rho for the area and multiplicty methods using kt jets
///

class DetermineTowerRho : public SubsysReco
{
 public:
  DetermineTowerRho(const std::string &name = "DetermineTowerRho");
  ~DetermineTowerRho() override;

  // standard Fun4All methods
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  // add rho method (Area or Multiplicity)
  void add_method(TowerRho::Method rho_method, std::string output = "")
  {
    // get method name
    std::string method_name = TowerRhov1::get_method_string(rho_method);

    // check if method already exists
    for (auto &method : _rho_methods)
    {
      if (method == rho_method)
      {
        std::cout << "WARNING: rho method " << method_name << " already exists" << std::endl;
        return;
      }
    }

    _rho_methods.push_back(rho_method);

    // if no output name is specified, use default
    if (output == "")
    {
      output = "TowerRho_" + method_name;
    }

    // add output name
    _outputs.push_back(output);
  }

  // inputs for background jets
  void add_tower_input(JetInput *input) { _inputs.push_back(input); }

  void set_algo(Jet::ALGO algo) { m_bkgd_jet_algo = algo; }
  Jet::ALGO get_algo() { return m_bkgd_jet_algo; }

  // set the jet algorithm parameter used to calculate the background jets
  void set_par(float par) { m_par = par; }  // default is 0.4
  float get_par() { return m_par; }

  // set the absolute eta range for the background jet calculation // default is 1.1
  void set_tower_abs_eta(float abseta) { m_abs_tower_eta_range = abseta; }
  float get_tower_abs_eta() const { return m_abs_tower_eta_range; }

  // set the absolute eta range for the background jet calculation // default is 1.1 - m_par (set in the init_algo method)
  void set_jet_abs_eta(float abseta) { m_abs_jet_eta_range = abseta; }
  float get_jet_abs_eta() const { return m_abs_jet_eta_range; }

  // set the number of hardest towers to omit from the background jet calculation // default is 2
  void set_omit_nhardest(int omit_nhardest) { m_omit_nhardest = omit_nhardest; }
  int get_omit_nhardest() const { return m_omit_nhardest; }

  // set the minimum pT for towers and jets used in the background jet calculation
  void set_tower_min_pT(float min_pT) { m_tower_min_pT = min_pT; }  // default is 0.0
  float get_tower_min_pT() const { return m_tower_min_pT; }

  // set the minimum pT for towers and jets used in the background jet calculation
  void set_jet_min_pT(float min_pT) { m_jet_min_pT = min_pT; }  // default is 1.0
  float get_jet_min_pT() const { return m_jet_min_pT; }

  // set the ghost area for the background jet calculation
  void set_ghost_area(float ghost_area) { m_ghost_area = ghost_area; }  // default is 0.01
  float get_ghost_area() const { return m_ghost_area; }

  // print settings
  void print_settings(std::ostream &os = std::cout) const;

  // tower cut
  // void do_tower_cut(bool b = true) { m_do_tower_cut = b; }

  // void set_tower_threshold(float threshold)
  // {
  // m_tower_threshold = threshold;
  // }

 private:
  // internal methods
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode, unsigned int ipos, const float rho, const float sigma);

  // variables
  std::vector<JetInput *> _inputs;
  std::vector<std::string> _outputs;
  std::vector<TowerRho::Method> _rho_methods;
  Jet::ALGO _algo{Jet::ALGO::KT};

  // int m_verbosity { 0 };

  // settings for the background jet calculation
  Jet::ALGO m_bkgd_jet_algo{Jet::ALGO::KT};  // default is KT

  float m_par{0.4};  // default is 0.4

  // eta ranges for the background jet calculation
  float m_abs_tower_eta_range{1.1};  // default is 1.1
  float m_abs_jet_eta_range{5.0};    // default is 1.1 - m_par (set in the init_algo method)

  float m_abs_ghost_eta{5.0};  // set in the init_algo method

  int m_omit_nhardest{2};

  float m_tower_min_pT{0.0};
  float m_jet_min_pT{0.0};

  float m_ghost_area{0.01};

  // tower threshold
  // bool m_do_tower_cut { false };
  // float m_tower_threshold { 0.0 };

  fastjet::JetAlgorithm get_fj_algo();
  std::string get_fj_algo_name();
  std::string m_fj_algo_name;
  std::vector<fastjet::PseudoJet> get_pseudojets(PHCompositeNode *topNode);

  float _percentile(const std::vector<float> &sorted_vec,
                    const float percentile,
                    const float nempty) const;
  void _median_stddev(const std::vector<float> &vec,
                      float n_empty_jets,
                      float &median,
                      float &std_dev) const;
};

#endif  // RHOBASE_DETERMINETOWERRHO_H

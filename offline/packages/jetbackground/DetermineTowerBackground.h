#ifndef JETBACKGROUND_DETERMINETOWERBACKGROUND_H
#define JETBACKGROUND_DETERMINETOWERBACKGROUND_H

//===========================================================
/// \file DetermineTowerBackground.h
/// \brief UE background calculator
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

// system includes
#include <jetbase/Jet.h>
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;

/// \class DetermineTowerBackground
///
/// \brief UE background calculator
///
/// This module constructs dE/deta vs. eta and v2 estimates given an
/// (unsubtracted) set of calorimeter towers and possible a set of
/// exclusion jets (seeds)
///
class DetermineTowerBackground : public SubsysReco
{
 public:
  DetermineTowerBackground(const std::string &name = "DetermineTowerBackground");
  ~DetermineTowerBackground() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void SetBackgroundOutputName(const std::string &name) { _backgroundName = name; }
  void SetSeedType(int seed_type) { _seed_type = seed_type; }
  void SetFlow(int do_flow) { _do_flow = do_flow; };

  void SetSeedJetD(float D) { _seed_jet_D = D; };
  void SetSeedJetPt(float pt) { _seed_jet_pt = pt; };
  void SetSeedMaxConst(float max_const) { _seed_max_const = max_const; };

  void UseReweighting(bool do_reweight ) {  _do_reweight = do_reweight; }

  void set_towerinfo(bool use_towerinfo)
  {
    m_use_towerinfo = use_towerinfo;
  }
  void set_towerNodePrefix(const std::string &prefix)
  {
    m_towerNodePrefix = prefix;
    return;
  }

 private:

  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int _do_flow{0};
  float _v2{0};
  float _Psi2{0};
  std::vector<std::vector<float> > _UE;
  int _nStrips{0};
  int _nTowers{0};

  int _HCAL_NETA{-1};
  int _HCAL_NPHI{-1};

  
  std::vector<std::vector<float> > _EMCAL_E;
  std::vector<std::vector<float> > _IHCAL_E;
  std::vector<std::vector<float> > _OHCAL_E;

  std::vector<std::vector<int> > _EMCAL_ISBAD;
  std::vector<std::vector<int> > _IHCAL_ISBAD;
  std::vector<std::vector<int> > _OHCAL_ISBAD;

  // 1-D energies vs. phi (integrated over eta strips with complete
  // phi coverage, and all layers)
  std::vector<float> _FULLCALOFLOW_PHI_E;
  std::vector<float> _FULLCALOFLOW_PHI_VAL;

  bool _do_reweight{true}; // flag to indicate if reweighting is used
  std::vector<float> _EMCAL_PHI_WEIGHTS;
  std::vector<float> _IHCAL_PHI_WEIGHTS;
  std::vector<float> _OHCAL_PHI_WEIGHTS;

  std::string _backgroundName{"TestTowerBackground"};

  int _seed_type{0};
  float _seed_jet_D{4.0};
  float _seed_max_const{3.0};
  float _seed_jet_pt{7.0};

  std::vector<float> _seed_eta;
  std::vector<float> _seed_phi;

  Jet::PROPERTY _index_SeedD{};
  Jet::PROPERTY _index_SeedItr{};

  bool m_use_towerinfo{false};
  bool _is_flow_failure{false};
  bool _reweight_failed{false};

  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;
};

#endif

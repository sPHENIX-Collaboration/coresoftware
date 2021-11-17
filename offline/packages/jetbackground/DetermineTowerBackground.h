#ifndef JETBACKGROUND_DETERMINETOWERBACKGROUND_H
#define JETBACKGROUND_DETERMINETOWERBACKGROUND_H

//===========================================================
/// \file DetermineTowerBackground.h
/// \brief UE background calculator
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

// system includes
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

 private:
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int _do_flow;
  float _v2;
  float _Psi2;
  std::vector<std::vector<float> > _UE;
  int _nStrips;
  int _nTowers;

  int _HCAL_NETA;
  int _HCAL_NPHI;

  std::vector<std::vector<float> > _EMCAL_E;
  std::vector<std::vector<float> > _IHCAL_E;
  std::vector<std::vector<float> > _OHCAL_E;

  // 1-D energies vs. phi (integrated over eta strips with complete
  // phi coverage, and all layers)
  std::vector<float> _FULLCALOFLOW_PHI_E;
  std::vector<float> _FULLCALOFLOW_PHI_VAL;

  std::string _backgroundName;

  int _seed_type;
  float _seed_jet_D;
  float _seed_jet_pt;

  std::vector<float> _seed_eta;
  std::vector<float> _seed_phi;
};

#endif

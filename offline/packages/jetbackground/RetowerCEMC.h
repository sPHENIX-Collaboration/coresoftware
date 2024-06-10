#ifndef JETBACKGROUND_RETOWERCEMC_H
#define JETBACKGROUND_RETOWERCEMC_H

//===========================================================
/// \file RetowerCEMC.h
/// \brief creates 0.1x0.1 towerized CEMC container
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

// system includes
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;

/// \class RetowerCEMC
///
/// \brief creates 0.1x0.1-towerized CEMC container
///
/// Using the existing CEMC tower collection, creates a new
/// 0.1x0.1-towerized version appropriate for HI jet reconstruction
///
class RetowerCEMC : public SubsysReco
{
 public:
  RetowerCEMC(const std::string &name = "RetowerCEMC");
  ~RetowerCEMC() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void SetEnergyDistribution(int val) { _WEIGHTED_ENERGY_DISTRIBUTION = val; }
  void set_towerinfo(bool use_towerinfo)
  {
    m_use_towerinfo = use_towerinfo;
  }
  void set_frac_cut(float frac_cut)
  {
    _FRAC_CUT = frac_cut;
  }
  void set_towerNodePrefix(const std::string &prefix)
  {
    m_towerNodePrefix = prefix;
    return;
  }

 private:
  int CreateNode(PHCompositeNode *topNode);
  int _WEIGHTED_ENERGY_DISTRIBUTION{1};
  int _NETA{-1};
  int _NPHI{-1};
  float _FRAC_CUT{0.5};
  bool m_use_towerinfo{false};
  std::vector<std::vector<float> > _EMCAL_RETOWER_E;
  std::vector<std::vector<int> > _EMCAL_RETOWER_T;
  std::vector<std::vector<float> > _EMCAL_RETOWER_MASKED_A;
  std::vector<std::vector<float> > _EMCAL_RETOWER_TOTAL_A;
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string EMTowerName;
  std::string IHTowerName;
  std::string EMRetowerName;
};

#endif

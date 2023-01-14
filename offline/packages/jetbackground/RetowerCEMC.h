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

 private:
  int CreateNode(PHCompositeNode *topNode);
  int _WEIGHTED_ENERGY_DISTRIBUTION;
  int _NETA;
  int _NPHI;
  std::vector<std::vector<float> > _EMCAL_RETOWER_E;
};

#endif

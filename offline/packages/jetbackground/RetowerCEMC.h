#ifndef JETBACKGROUND_RETOWERCEMC_H
#define JETBACKGROUND_RETOWERCEMC_H

//===========================================================
/// \file RetowerCEMC.h
/// \brief creates 0.1x0.1 towerized CEMC container
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

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
  virtual ~RetowerCEMC();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  int CreateNode(PHCompositeNode *topNode);

  int _NETA;
  int _NPHI;
  std::vector< std::vector<float> > _EMCAL_RETOWER_E;

};

#endif

#ifndef __CALOTRIGGERSIM_H__
#define __CALOTRIGGERSIM_H__

//===========================================================
/// \file CaloTriggerSim.h
/// \brief simple trigger emulation
/// \author Dennis V. Perepelitsa
//===========================================================

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

// standard includes
#include <vector>

// forward declarations
class PHCompositeNode;

/// \class CaloTriggerSim
///
/// \brief simple trigger emulation
///
/// This module constructs calo-based trigger signatures
///
class CaloTriggerSim : public SubsysReco
{
 public:
  CaloTriggerSim(const std::string &name = "CaloTriggerSim");
  virtual ~CaloTriggerSim();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int _EMCAL_1x1_NETA;
  int _EMCAL_1x1_NPHI;
  std::vector<std::vector<float> > _EMCAL_1x1_MAP;

  int _EMCAL_2x2_NETA;
  int _EMCAL_2x2_NPHI;
  std::vector<std::vector<float> > _EMCAL_2x2_MAP;

  int _EMCAL_4x4_NETA;
  int _EMCAL_4x4_NPHI;
  std::vector<std::vector<float> > _EMCAL_4x4_MAP;

  float _EMCAL_2x2_BEST_E;
  float _EMCAL_2x2_BEST_PHI;
  float _EMCAL_2x2_BEST_ETA;

  float _EMCAL_4x4_BEST_E;
  float _EMCAL_4x4_BEST_PHI;
  float _EMCAL_4x4_BEST_ETA;

  //int verbosity;
};

#endif  // __CALOTRIGGERSIM_H__

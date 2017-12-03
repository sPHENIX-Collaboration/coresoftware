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

  void set_truncation( int emulate_truncation );
  float truncate_8bit( float raw_E );

 private:
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int _emulate_truncation;

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

  float _EMCAL_4x4_BEST2_E;
  float _EMCAL_4x4_BEST2_PHI;
  float _EMCAL_4x4_BEST2_ETA;

  // needed since phi ranges are potentially different in EMCal vs. HCal
  float _FULLCALO_PHI_START;
  float _FULLCALO_PHI_END;

  // full calo (based on 0.1x0.1 HCal towers) limits and maps
  int _FULLCALO_0p1x0p1_NETA;
  int _FULLCALO_0p1x0p1_NPHI;
  std::vector<std::vector<float> > _FULLCALO_0p1x0p1_MAP;

  int _FULLCALO_0p2x0p2_NETA;
  int _FULLCALO_0p2x0p2_NPHI;
  std::vector<std::vector<float> > _FULLCALO_0p2x0p2_MAP;

  int _FULLCALO_0p4x0p4_NETA;
  int _FULLCALO_0p4x0p4_NPHI;
  std::vector<std::vector<float> > _FULLCALO_0p4x0p4_MAP;

  int _FULLCALO_0p6x0p6_NETA;
  int _FULLCALO_0p6x0p6_NPHI;
  std::vector<std::vector<float> > _FULLCALO_0p6x0p6_MAP;

  int _FULLCALO_0p8x0p8_NETA;
  int _FULLCALO_0p8x0p8_NPHI;
  std::vector<std::vector<float> > _FULLCALO_0p8x0p8_MAP;

  int _FULLCALO_1p0x1p0_NETA;
  int _FULLCALO_1p0x1p0_NPHI;
  std::vector<std::vector<float> > _FULLCALO_1p0x1p0_MAP;
  
  // highest full calo window energies
  float _FULLCALO_0p2x0p2_BEST_E;
  float _FULLCALO_0p2x0p2_BEST_PHI;
  float _FULLCALO_0p2x0p2_BEST_ETA;

  float _FULLCALO_0p4x0p4_BEST_E;
  float _FULLCALO_0p4x0p4_BEST_PHI;
  float _FULLCALO_0p4x0p4_BEST_ETA;

  float _FULLCALO_0p6x0p6_BEST_E;
  float _FULLCALO_0p6x0p6_BEST_PHI;
  float _FULLCALO_0p6x0p6_BEST_ETA;

  float _FULLCALO_0p8x0p8_BEST_E;
  float _FULLCALO_0p8x0p8_BEST_PHI;
  float _FULLCALO_0p8x0p8_BEST_ETA;

  float _FULLCALO_1p0x1p0_BEST_E;
  float _FULLCALO_1p0x1p0_BEST_PHI;
  float _FULLCALO_1p0x1p0_BEST_ETA;

  //int verbosity;
};

#endif  // __CALOTRIGGERSIM_H__

#ifndef TRIGGER_CALOTRIGGERSIM_H
#define TRIGGER_CALOTRIGGERSIM_H

//===========================================================
/// \file CaloTriggerSim.h
/// \brief simple trigger emulation
/// \author Dennis V. Perepelitsa
//===========================================================

// sPHENIX includes
#include <fun4all/SubsysReco.h>

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
  virtual ~CaloTriggerSim() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Reset();

  void set_truncation(const int emulate_truncation) {m_EmulateTruncationFlag = emulate_truncation;}
  double truncate_8bit(const double raw_E) const;

 private:
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int m_EmulateTruncationFlag;

  int m_EMCAL_1x1_NETA;
  int m_EMCAL_1x1_NPHI;

  int m_EMCAL_2x2_NETA;
  int m_EMCAL_2x2_NPHI;

  int m_EMCAL_4x4_NETA;
  int m_EMCAL_4x4_NPHI;

  float m_EMCAL_2x2_BEST_E;
  float m_EMCAL_2x2_BEST_PHI;
  float m_EMCAL_2x2_BEST_ETA;

  float m_EMCAL_4x4_BEST_E;
  float m_EMCAL_4x4_BEST_PHI;
  float m_EMCAL_4x4_BEST_ETA;

  float m_EMCAL_4x4_BEST2_E;
  float m_EMCAL_4x4_BEST2_PHI;
  float m_EMCAL_4x4_BEST2_ETA;

  // needed since phi ranges are potentially different in EMCal vs. HCal
  float m_FULLCALO_PHI_START;
  float m_FULLCALO_PHI_END;

  // full calo (based on 0.1x0.1 HCal towers) limits and maps
  int m_FULLCALO_0p1x0p1_NETA;
  int m_FULLCALO_0p1x0p1_NPHI;

  int m_FULLCALO_0p2x0p2_NETA;
  int m_FULLCALO_0p2x0p2_NPHI;

  int m_FULLCALO_0p4x0p4_NETA;
  int m_FULLCALO_0p4x0p4_NPHI;

  int m_FULLCALO_0p6x0p6_NETA;
  int m_FULLCALO_0p6x0p6_NPHI;

  int m_FULLCALO_0p8x0p8_NETA;
  int m_FULLCALO_0p8x0p8_NPHI;

  int m_FULLCALO_1p0x1p0_NETA;
  int m_FULLCALO_1p0x1p0_NPHI;

  // highest full calo window energies
  float m_FULLCALO_0p2x0p2_BEST_E;
  float m_FULLCALO_0p2x0p2_BEST_PHI;
  float m_FULLCALO_0p2x0p2_BEST_ETA;

  float m_FULLCALO_0p4x0p4_BEST_E;
  float m_FULLCALO_0p4x0p4_BEST_PHI;
  float m_FULLCALO_0p4x0p4_BEST_ETA;

  float m_FULLCALO_0p6x0p6_BEST_E;
  float m_FULLCALO_0p6x0p6_BEST_PHI;
  float m_FULLCALO_0p6x0p6_BEST_ETA;

  float m_FULLCALO_0p8x0p8_BEST_E;
  float m_FULLCALO_0p8x0p8_BEST_PHI;
  float m_FULLCALO_0p8x0p8_BEST_ETA;

  float m_FULLCALO_1p0x1p0_BEST_E;
  float m_FULLCALO_1p0x1p0_BEST_PHI;
  float m_FULLCALO_1p0x1p0_BEST_ETA;

  std::vector<std::vector<double> > m_EMCAL_1x1_MAP;
  std::vector<std::vector<double> > m_EMCAL_2x2_MAP;
  std::vector<std::vector<double> > m_EMCAL_4x4_MAP;
  std::vector<std::vector<double> > m_FULLCALO_0p1x0p1_MAP;
  std::vector<std::vector<double> > m_FULLCALO_0p2x0p2_MAP;
  std::vector<std::vector<double> > m_FULLCALO_0p4x0p4_MAP;
  std::vector<std::vector<double> > m_FULLCALO_0p6x0p6_MAP;
  std::vector<std::vector<double> > m_FULLCALO_0p8x0p8_MAP;
  std::vector<std::vector<double> > m_FULLCALO_1p0x1p0_MAP;
};

#endif  // TRIGGER_CALOTRIGGERSIM_H

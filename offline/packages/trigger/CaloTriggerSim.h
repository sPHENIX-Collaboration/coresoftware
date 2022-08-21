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
#include <cmath>
#include <string>
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
  ~CaloTriggerSim() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_truncation(const int emulate_truncation) { m_EmulateTruncationFlag = emulate_truncation; }
  double truncate_8bit(const double raw_E) const;

 private:
  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  int m_EmulateTruncationFlag = 0;

  int m_EMCAL_1x1_NETA = -1;
  int m_EMCAL_1x1_NPHI = -1;

  int m_EMCAL_2x2_NETA = -1;
  int m_EMCAL_2x2_NPHI = -1;

  int m_EMCAL_4x4_NETA = -1;
  int m_EMCAL_4x4_NPHI = -1;

  float m_EMCAL_2x2_BEST_E = 0.;
  float m_EMCAL_2x2_BEST_PHI = 0.;
  float m_EMCAL_2x2_BEST_ETA = 0.;

  float m_EMCAL_4x4_BEST_E = 0.;
  float m_EMCAL_4x4_BEST_PHI = 0.;
  float m_EMCAL_4x4_BEST_ETA = 0.;

  float m_EMCAL_4x4_BEST2_E = 0.;
  float m_EMCAL_4x4_BEST2_PHI = 0.;
  float m_EMCAL_4x4_BEST2_ETA = 0.;

  // needed since phi ranges are potentially different in EMCal vs. HCal
  float m_FULLCALO_PHI_START = 0.;
  float m_FULLCALO_PHI_END = 2 * M_PI;

  // full calo (based on 0.1x0.1 HCal towers) limits and maps
  int m_FULLCALO_0p1x0p1_NETA = -1;
  int m_FULLCALO_0p1x0p1_NPHI = -1;

  int m_FULLCALO_0p2x0p2_NETA = -1;
  int m_FULLCALO_0p2x0p2_NPHI = -1;

  int m_FULLCALO_0p4x0p4_NETA = -1;
  int m_FULLCALO_0p4x0p4_NPHI = -1;

  int m_FULLCALO_0p6x0p6_NETA = -1;
  int m_FULLCALO_0p6x0p6_NPHI = -1;

  int m_FULLCALO_0p8x0p8_NETA = -1;
  int m_FULLCALO_0p8x0p8_NPHI = -1;

  int m_FULLCALO_1p0x1p0_NETA = -1;
  int m_FULLCALO_1p0x1p0_NPHI = -1;

  // highest full calo window energies
  float m_FULLCALO_0p2x0p2_BEST_E = 0.;
  float m_FULLCALO_0p2x0p2_BEST_PHI = 0.;
  float m_FULLCALO_0p2x0p2_BEST_ETA = 0.;

  float m_FULLCALO_0p4x0p4_BEST_E = 0.;
  float m_FULLCALO_0p4x0p4_BEST_PHI = 0.;
  float m_FULLCALO_0p4x0p4_BEST_ETA = 0.;

  float m_FULLCALO_0p6x0p6_BEST_E = 0.;
  float m_FULLCALO_0p6x0p6_BEST_PHI = 0.;
  float m_FULLCALO_0p6x0p6_BEST_ETA = 0.;

  float m_FULLCALO_0p8x0p8_BEST_E = 0.;
  float m_FULLCALO_0p8x0p8_BEST_PHI = 0.;
  float m_FULLCALO_0p8x0p8_BEST_ETA = 0.;

  float m_FULLCALO_1p0x1p0_BEST_E = 0.;
  float m_FULLCALO_1p0x1p0_BEST_PHI = 0.;
  float m_FULLCALO_1p0x1p0_BEST_ETA = 0.;

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

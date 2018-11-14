#ifndef TRIGGER_CALOTRIGGERINFOV1_H
#define TRIGGER_CALOTRIGGERINFOV1_H

#include "CaloTriggerInfo.h"

#include <phool/PHObject.h>

class CaloTriggerInfo_v1 : public CaloTriggerInfo
{
 public:
  CaloTriggerInfo_v1();
  virtual ~CaloTriggerInfo_v1();

  void identify(std::ostream &os = std::cout) const;
  void Reset() {}
  int isValid() const { return 1; }

  // EMCal 2x2
  void set_best_EMCal_2x2_E(const float E) { m_EMCAL_2x2_BEST_E = E; }
  void set_best_EMCal_2x2_eta(const float eta) { m_EMCAL_2x2_BEST_ETA = eta; }
  void set_best_EMCal_2x2_phi(const float phi) { m_EMCAL_2x2_BEST_PHI = phi; }

  float get_best_EMCal_2x2_E() const { return m_EMCAL_2x2_BEST_E; }
  float get_best_EMCal_2x2_eta() const { return m_EMCAL_2x2_BEST_ETA; }
  float get_best_EMCal_2x2_phi() const { return m_EMCAL_2x2_BEST_PHI; }

  // EMCal 4x4
  void set_best_EMCal_4x4_E(const float E) { m_EMCAL_4x4_BEST_E = E; }
  void set_best_EMCal_4x4_eta(const float eta) { m_EMCAL_4x4_BEST_ETA = eta; }
  void set_best_EMCal_4x4_phi(const float phi) { m_EMCAL_4x4_BEST_PHI = phi; }

  float get_best_EMCal_4x4_E() const { return m_EMCAL_4x4_BEST_E; }
  float get_best_EMCal_4x4_eta() const { return m_EMCAL_4x4_BEST_ETA; }
  float get_best_EMCal_4x4_phi() const { return m_EMCAL_4x4_BEST_PHI; }

  // 2nd best EMCal 4x4
  void set_best2_EMCal_4x4_E(const float E) { m_EMCAL_4x4_BEST2_E = E; }
  void set_best2_EMCal_4x4_eta(const float eta) { m_EMCAL_4x4_BEST2_ETA = eta; }
  void set_best2_EMCal_4x4_phi(const float phi) { m_EMCAL_4x4_BEST2_PHI = phi; }

  float get_best2_EMCal_4x4_E() const { return m_EMCAL_4x4_BEST2_E; }
  float get_best2_EMCal_4x4_eta() const { return m_EMCAL_4x4_BEST2_ETA; }
  float get_best2_EMCal_4x4_phi() const { return m_EMCAL_4x4_BEST2_PHI; }

  // FullCalo 0.2x0.2
  void set_best_FullCalo_0p2x0p2_E(const float E) { m_FULLCALO_0p2x0p2_BEST_E = E; }
  void set_best_FullCalo_0p2x0p2_eta(const float eta) { m_FULLCALO_0p2x0p2_BEST_ETA = eta; }
  void set_best_FullCalo_0p2x0p2_phi(const float phi) { m_FULLCALO_0p2x0p2_BEST_PHI = phi; }

  float get_best_FullCalo_0p2x0p2_E() const { return m_FULLCALO_0p2x0p2_BEST_E; }
  float get_best_FullCalo_0p2x0p2_eta() const { return m_FULLCALO_0p2x0p2_BEST_ETA; }
  float get_best_FullCalo_0p2x0p2_phi() const { return m_FULLCALO_0p2x0p2_BEST_PHI; }

  // FullCalo 0.4x0.4
  void set_best_FullCalo_0p4x0p4_E(const float E) { m_FULLCALO_0p4x0p4_BEST_E = E; }
  void set_best_FullCalo_0p4x0p4_eta(const float eta) { m_FULLCALO_0p4x0p4_BEST_ETA = eta; }
  void set_best_FullCalo_0p4x0p4_phi(const float phi) { m_FULLCALO_0p4x0p4_BEST_PHI = phi; }

  float get_best_FullCalo_0p4x0p4_E() const { return m_FULLCALO_0p4x0p4_BEST_E; }
  float get_best_FullCalo_0p4x0p4_eta() const { return m_FULLCALO_0p4x0p4_BEST_ETA; }
  float get_best_FullCalo_0p4x0p4_phi() const { return m_FULLCALO_0p4x0p4_BEST_PHI; }

  // FullCalo 0.6x0.6
  void set_best_FullCalo_0p6x0p6_E(const float E) { m_FULLCALO_0p6x0p6_BEST_E = E; }
  void set_best_FullCalo_0p6x0p6_eta(const float eta) { m_FULLCALO_0p6x0p6_BEST_ETA = eta; }
  void set_best_FullCalo_0p6x0p6_phi(const float phi) { m_FULLCALO_0p6x0p6_BEST_PHI = phi; }

  float get_best_FullCalo_0p6x0p6_E() const { return m_FULLCALO_0p6x0p6_BEST_E; }
  float get_best_FullCalo_0p6x0p6_eta() const { return m_FULLCALO_0p6x0p6_BEST_ETA; }
  float get_best_FullCalo_0p6x0p6_phi() const { return m_FULLCALO_0p6x0p6_BEST_PHI; }

  // FullCalo 0.8x0.8
  void set_best_FullCalo_0p8x0p8_E(const float E) { m_FULLCALO_0p8x0p8_BEST_E = E; }
  void set_best_FullCalo_0p8x0p8_eta(const float eta) { m_FULLCALO_0p8x0p8_BEST_ETA = eta; }
  void set_best_FullCalo_0p8x0p8_phi(const float phi) { m_FULLCALO_0p8x0p8_BEST_PHI = phi; }

  float get_best_FullCalo_0p8x0p8_E() const { return m_FULLCALO_0p8x0p8_BEST_E; }
  float get_best_FullCalo_0p8x0p8_eta() const { return m_FULLCALO_0p8x0p8_BEST_ETA; }
  float get_best_FullCalo_0p8x0p8_phi() const { return m_FULLCALO_0p8x0p8_BEST_PHI; }

  // FullCalo 1.0x1.0
  void set_best_FullCalo_1p0x1p0_E(const float E) { m_FULLCALO_1p0x1p0_BEST_E = E; }
  void set_best_FullCalo_1p0x1p0_eta(const float eta) { m_FULLCALO_1p0x1p0_BEST_ETA = eta; }
  void set_best_FullCalo_1p0x1p0_phi(const float phi) { m_FULLCALO_1p0x1p0_BEST_PHI = phi; }

  float get_best_FullCalo_1p0x1p0_E() const { return m_FULLCALO_1p0x1p0_BEST_E; }
  float get_best_FullCalo_1p0x1p0_eta() const { return m_FULLCALO_1p0x1p0_BEST_ETA; }
  float get_best_FullCalo_1p0x1p0_phi() const { return m_FULLCALO_1p0x1p0_BEST_PHI; }

 private:
  float m_EMCAL_2x2_BEST_E;
  float m_EMCAL_2x2_BEST_ETA;
  float m_EMCAL_2x2_BEST_PHI;

  float m_EMCAL_4x4_BEST_E;
  float m_EMCAL_4x4_BEST_ETA;
  float m_EMCAL_4x4_BEST_PHI;

  float m_EMCAL_4x4_BEST2_E;
  float m_EMCAL_4x4_BEST2_ETA;
  float m_EMCAL_4x4_BEST2_PHI;

  float m_FULLCALO_0p2x0p2_BEST_E;
  float m_FULLCALO_0p2x0p2_BEST_ETA;
  float m_FULLCALO_0p2x0p2_BEST_PHI;

  float m_FULLCALO_0p4x0p4_BEST_E;
  float m_FULLCALO_0p4x0p4_BEST_ETA;
  float m_FULLCALO_0p4x0p4_BEST_PHI;

  float m_FULLCALO_0p6x0p6_BEST_E;
  float m_FULLCALO_0p6x0p6_BEST_ETA;
  float m_FULLCALO_0p6x0p6_BEST_PHI;

  float m_FULLCALO_0p8x0p8_BEST_E;
  float m_FULLCALO_0p8x0p8_BEST_ETA;
  float m_FULLCALO_0p8x0p8_BEST_PHI;

  float m_FULLCALO_1p0x1p0_BEST_E;
  float m_FULLCALO_1p0x1p0_BEST_ETA;
  float m_FULLCALO_1p0x1p0_BEST_PHI;

  ClassDef(CaloTriggerInfo_v1, 3);
};

#endif  // __CALOTRIGGERINFO_V1_H__

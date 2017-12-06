#ifndef __CALOTRIGGERINFO_V1_H__
#define __CALOTRIGGERINFO_V1_H__

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
  void set_best_EMCal_2x2_E(float E) { _EMCAL_2x2_BEST_E = E; }
  void set_best_EMCal_2x2_eta(float eta) { _EMCAL_2x2_BEST_ETA = eta; }
  void set_best_EMCal_2x2_phi(float phi) { _EMCAL_2x2_BEST_PHI = phi; }

  float get_best_EMCal_2x2_E() { return _EMCAL_2x2_BEST_E; }
  float get_best_EMCal_2x2_eta() { return _EMCAL_2x2_BEST_ETA; }
  float get_best_EMCal_2x2_phi() { return _EMCAL_2x2_BEST_PHI; }

  // EMCal 4x4
  void set_best_EMCal_4x4_E(float E) { _EMCAL_4x4_BEST_E = E; }
  void set_best_EMCal_4x4_eta(float eta) { _EMCAL_4x4_BEST_ETA = eta; }
  void set_best_EMCal_4x4_phi(float phi) { _EMCAL_4x4_BEST_PHI = phi; }

  float get_best_EMCal_4x4_E() { return _EMCAL_4x4_BEST_E; }
  float get_best_EMCal_4x4_eta() { return _EMCAL_4x4_BEST_ETA; }
  float get_best_EMCal_4x4_phi() { return _EMCAL_4x4_BEST_PHI; }

  // 2nd best EMCal 4x4
  void set_best2_EMCal_4x4_E(float E) { _EMCAL_4x4_BEST2_E = E; }
  void set_best2_EMCal_4x4_eta(float eta) { _EMCAL_4x4_BEST2_ETA = eta; }
  void set_best2_EMCal_4x4_phi(float phi) { _EMCAL_4x4_BEST2_PHI = phi; }

  float get_best2_EMCal_4x4_E() { return _EMCAL_4x4_BEST2_E; }
  float get_best2_EMCal_4x4_eta() { return _EMCAL_4x4_BEST2_ETA; }
  float get_best2_EMCal_4x4_phi() { return _EMCAL_4x4_BEST2_PHI; }

  // FullCalo 0.2x0.2
  void set_best_FullCalo_0p2x0p2_E(float E) { _FULLCALO_0p2x0p2_BEST_E = E; }
  void set_best_FullCalo_0p2x0p2_eta(float eta) { _FULLCALO_0p2x0p2_BEST_ETA = eta; }
  void set_best_FullCalo_0p2x0p2_phi(float phi) { _FULLCALO_0p2x0p2_BEST_PHI = phi; }

  float get_best_FullCalo_0p2x0p2_E() { return _FULLCALO_0p2x0p2_BEST_E; }
  float get_best_FullCalo_0p2x0p2_eta() { return _FULLCALO_0p2x0p2_BEST_ETA; }
  float get_best_FullCalo_0p2x0p2_phi() { return _FULLCALO_0p2x0p2_BEST_PHI; }

  // FullCalo 0.4x0.4
  void set_best_FullCalo_0p4x0p4_E(float E) { _FULLCALO_0p4x0p4_BEST_E = E; }
  void set_best_FullCalo_0p4x0p4_eta(float eta) { _FULLCALO_0p4x0p4_BEST_ETA = eta; }
  void set_best_FullCalo_0p4x0p4_phi(float phi) { _FULLCALO_0p4x0p4_BEST_PHI = phi; }

  float get_best_FullCalo_0p4x0p4_E() { return _FULLCALO_0p4x0p4_BEST_E; }
  float get_best_FullCalo_0p4x0p4_eta() { return _FULLCALO_0p4x0p4_BEST_ETA; }
  float get_best_FullCalo_0p4x0p4_phi() { return _FULLCALO_0p4x0p4_BEST_PHI; }

  // FullCalo 0.6x0.6
  void set_best_FullCalo_0p6x0p6_E(float E) { _FULLCALO_0p6x0p6_BEST_E = E; }
  void set_best_FullCalo_0p6x0p6_eta(float eta) { _FULLCALO_0p6x0p6_BEST_ETA = eta; }
  void set_best_FullCalo_0p6x0p6_phi(float phi) { _FULLCALO_0p6x0p6_BEST_PHI = phi; }

  float get_best_FullCalo_0p6x0p6_E() { return _FULLCALO_0p6x0p6_BEST_E; }
  float get_best_FullCalo_0p6x0p6_eta() { return _FULLCALO_0p6x0p6_BEST_ETA; }
  float get_best_FullCalo_0p6x0p6_phi() { return _FULLCALO_0p6x0p6_BEST_PHI; }

  // FullCalo 0.8x0.8
  void set_best_FullCalo_0p8x0p8_E(float E) { _FULLCALO_0p8x0p8_BEST_E = E; }
  void set_best_FullCalo_0p8x0p8_eta(float eta) { _FULLCALO_0p8x0p8_BEST_ETA = eta; }
  void set_best_FullCalo_0p8x0p8_phi(float phi) { _FULLCALO_0p8x0p8_BEST_PHI = phi; }

  float get_best_FullCalo_0p8x0p8_E() { return _FULLCALO_0p8x0p8_BEST_E; }
  float get_best_FullCalo_0p8x0p8_eta() { return _FULLCALO_0p8x0p8_BEST_ETA; }
  float get_best_FullCalo_0p8x0p8_phi() { return _FULLCALO_0p8x0p8_BEST_PHI; }

  // FullCalo 1.0x1.0
  void set_best_FullCalo_1p0x1p0_E(float E) { _FULLCALO_1p0x1p0_BEST_E = E; }
  void set_best_FullCalo_1p0x1p0_eta(float eta) { _FULLCALO_1p0x1p0_BEST_ETA = eta; }
  void set_best_FullCalo_1p0x1p0_phi(float phi) { _FULLCALO_1p0x1p0_BEST_PHI = phi; }

  float get_best_FullCalo_1p0x1p0_E() { return _FULLCALO_1p0x1p0_BEST_E; }
  float get_best_FullCalo_1p0x1p0_eta() { return _FULLCALO_1p0x1p0_BEST_ETA; }
  float get_best_FullCalo_1p0x1p0_phi() { return _FULLCALO_1p0x1p0_BEST_PHI; }

 private:
  float _EMCAL_2x2_BEST_E;
  float _EMCAL_2x2_BEST_ETA;
  float _EMCAL_2x2_BEST_PHI;

  float _EMCAL_4x4_BEST_E;
  float _EMCAL_4x4_BEST_ETA;
  float _EMCAL_4x4_BEST_PHI;

  float _EMCAL_4x4_BEST2_E;
  float _EMCAL_4x4_BEST2_ETA;
  float _EMCAL_4x4_BEST2_PHI;

  float _FULLCALO_0p2x0p2_BEST_E;
  float _FULLCALO_0p2x0p2_BEST_ETA;
  float _FULLCALO_0p2x0p2_BEST_PHI;

  float _FULLCALO_0p4x0p4_BEST_E;
  float _FULLCALO_0p4x0p4_BEST_ETA;
  float _FULLCALO_0p4x0p4_BEST_PHI;

  float _FULLCALO_0p6x0p6_BEST_E;
  float _FULLCALO_0p6x0p6_BEST_ETA;
  float _FULLCALO_0p6x0p6_BEST_PHI;

  float _FULLCALO_0p8x0p8_BEST_E;
  float _FULLCALO_0p8x0p8_BEST_ETA;
  float _FULLCALO_0p8x0p8_BEST_PHI;

  float _FULLCALO_1p0x1p0_BEST_E;
  float _FULLCALO_1p0x1p0_BEST_ETA;
  float _FULLCALO_1p0x1p0_BEST_PHI;

  ClassDef(CaloTriggerInfo_v1, 2);
};

#endif  // __CALOTRIGGERINFO_V1_H__

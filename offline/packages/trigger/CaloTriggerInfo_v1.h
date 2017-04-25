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
  void set_best_2x2_E(float E) { _EMCAL_2x2_BEST_E = E; }
  void set_best_2x2_eta(float eta) { _EMCAL_2x2_BEST_ETA = eta; }
  void set_best_2x2_phi(float phi) { _EMCAL_2x2_BEST_PHI = phi; }
  float get_best_2x2_E() { return _EMCAL_2x2_BEST_E; }
  float get_best_2x2_eta() { return _EMCAL_2x2_BEST_ETA; }
  float get_best_2x2_phi() { return _EMCAL_2x2_BEST_PHI; }
  void set_best_4x4_E(float E) { _EMCAL_4x4_BEST_E = E; }
  void set_best_4x4_eta(float eta) { _EMCAL_4x4_BEST_ETA = eta; }
  void set_best_4x4_phi(float phi) { _EMCAL_4x4_BEST_PHI = phi; }
  float get_best_4x4_E() { return _EMCAL_4x4_BEST_E; }
  float get_best_4x4_eta() { return _EMCAL_4x4_BEST_ETA; }
  float get_best_4x4_phi() { return _EMCAL_4x4_BEST_PHI; }
 private:
  float _EMCAL_2x2_BEST_E;
  float _EMCAL_2x2_BEST_ETA;
  float _EMCAL_2x2_BEST_PHI;

  float _EMCAL_4x4_BEST_E;
  float _EMCAL_4x4_BEST_ETA;
  float _EMCAL_4x4_BEST_PHI;

  ClassDef(CaloTriggerInfo_v1, 1);
};

#endif  // __CALOTRIGGERINFO_V1_H__

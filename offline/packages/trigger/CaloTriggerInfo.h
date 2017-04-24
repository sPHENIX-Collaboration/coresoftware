#ifndef __CALOTRIGGERINFO_H__
#define __CALOTRIGGERINFO_H__

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class CaloTriggerInfo : public PHObject
{
 public:
  CaloTriggerInfo();
  virtual ~CaloTriggerInfo();

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

  ClassDef(CaloTriggerInfo, 1);
};

#endif  // __CALOTRIGGERINFO_H__

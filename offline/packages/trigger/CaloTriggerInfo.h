#ifndef __CALOTRIGGERINFO_H__
#define __CALOTRIGGERINFO_H__

#include <phool/PHObject.h>

class CaloTriggerInfo : public PHObject
{
 public:
  virtual ~CaloTriggerInfo() {};

  virtual void identify(std::ostream &os = std::cout) const { os << "CaloTriggerInfo base class" << std::endl; };
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  // EMCal 2x2
  virtual void set_best_EMCal_2x2_E(float E) {}
  virtual void set_best_EMCal_2x2_eta(float eta) {}
  virtual void set_best_EMCal_2x2_phi(float phi) {}

  virtual float get_best_EMCal_2x2_E() { return 0; }
  virtual float get_best_EMCal_2x2_eta() { return 0; }
  virtual float get_best_EMCal_2x2_phi() { return 0; }

  // EMCal 4x4
  virtual void set_best_EMCal_4x4_E(float E) {}
  virtual void set_best_EMCal_4x4_eta(float eta) {}
  virtual void set_best_EMCal_4x4_phi(float phi) {}

  virtual float get_best_EMCal_4x4_E() { return 0; }
  virtual float get_best_EMCal_4x4_eta() { return 0; }
  virtual float get_best_EMCal_4x4_phi() { return 0; }
  
  // 2nd best
  virtual void set_best2_EMCal_4x4_E(float E) {}
  virtual void set_best2_EMCal_4x4_eta(float eta) {}
  virtual void set_best2_EMCal_4x4_phi(float phi) {}

  virtual float get_best2_EMCal_4x4_E() { return 0; }
  virtual float get_best2_EMCal_4x4_eta() { return 0; }
  virtual float get_best2_EMCal_4x4_phi() { return 0; }

  // FullCalo 0.2x0.2
  virtual void set_best_FullCalo_0p2x0p2_E(float E) {}
  virtual void set_best_FullCalo_0p2x0p2_eta(float eta) {}
  virtual void set_best_FullCalo_0p2x0p2_phi(float phi) {}

  virtual float get_best_FullCalo_0p2x0p2_E() { return 0; }
  virtual float get_best_FullCalo_0p2x0p2_eta() { return 0; }
  virtual float get_best_FullCalo_0p2x0p2_phi() { return 0; }

  // FullCalo 0.4x0.4
  virtual void set_best_FullCalo_0p4x0p4_E(float E) {}
  virtual void set_best_FullCalo_0p4x0p4_eta(float eta) {}
  virtual void set_best_FullCalo_0p4x0p4_phi(float phi) {}

  virtual float get_best_FullCalo_0p4x0p4_E() { return 0; }
  virtual float get_best_FullCalo_0p4x0p4_eta() { return 0; }
  virtual float get_best_FullCalo_0p4x0p4_phi() { return 0; }

  // FullCalo 0.6x0.6
  virtual void set_best_FullCalo_0p6x0p6_E(float E) {}
  virtual void set_best_FullCalo_0p6x0p6_eta(float eta) {}
  virtual void set_best_FullCalo_0p6x0p6_phi(float phi) {}

  virtual float get_best_FullCalo_0p6x0p6_E() { return 0; }
  virtual float get_best_FullCalo_0p6x0p6_eta() { return 0; }
  virtual float get_best_FullCalo_0p6x0p6_phi() { return 0; }

  // FullCalo 0.8x0.8
  virtual void set_best_FullCalo_0p8x0p8_E(float E) {}
  virtual void set_best_FullCalo_0p8x0p8_eta(float eta) {}
  virtual void set_best_FullCalo_0p8x0p8_phi(float phi) {}

  virtual float get_best_FullCalo_0p8x0p8_E() { return 0; }
  virtual float get_best_FullCalo_0p8x0p8_eta() { return 0; }
  virtual float get_best_FullCalo_0p8x0p8_phi() { return 0; }

  // FullCalo 1.0x1.0
  virtual void set_best_FullCalo_1p0x1p0_E(float E) {}
  virtual void set_best_FullCalo_1p0x1p0_eta(float eta) {}
  virtual void set_best_FullCalo_1p0x1p0_phi(float phi) {}

  virtual float get_best_FullCalo_1p0x1p0_E() { return 0; }
  virtual float get_best_FullCalo_1p0x1p0_eta() { return 0; }
  virtual float get_best_FullCalo_1p0x1p0_phi() { return 0; }

 protected:
  CaloTriggerInfo() {}

 private:

  ClassDef(CaloTriggerInfo, 1);
};

#endif  // __CALOTRIGGERINFO_H__

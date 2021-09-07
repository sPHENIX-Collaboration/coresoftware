#ifndef TRIGGER_CALOTRIGGERINFO_H
#define TRIGGER_CALOTRIGGERINFO_H

#include <phool/PHObject.h>

class CaloTriggerInfo : public PHObject
{
 public:
  ~CaloTriggerInfo() override{};

  void identify(std::ostream &os = std::cout) const override { os << "CaloTriggerInfo base class" << std::endl; };
  int isValid() const override { return 0; }

  // EMCal 2x2
  virtual void set_best_EMCal_2x2_E(const float /*E*/) {}
  virtual void set_best_EMCal_2x2_eta(const float /*eta*/) {}
  virtual void set_best_EMCal_2x2_phi(const float /*phi*/) {}

  virtual float get_best_EMCal_2x2_E() const { return 0; }
  virtual float get_best_EMCal_2x2_eta() const { return 0; }
  virtual float get_best_EMCal_2x2_phi() const { return 0; }

  // EMCal 4x4
  virtual void set_best_EMCal_4x4_E(const float /*E*/) {}
  virtual void set_best_EMCal_4x4_eta(const float /*eta*/) {}
  virtual void set_best_EMCal_4x4_phi(const float /*phi*/) {}

  virtual float get_best_EMCal_4x4_E() const { return 0; }
  virtual float get_best_EMCal_4x4_eta() const { return 0; }
  virtual float get_best_EMCal_4x4_phi() const { return 0; }

  // 2nd best
  virtual void set_best2_EMCal_4x4_E(const float /*E*/) {}
  virtual void set_best2_EMCal_4x4_eta(const float /*eta*/) {}
  virtual void set_best2_EMCal_4x4_phi(const float /*phi*/) {}

  virtual float get_best2_EMCal_4x4_E() const { return 0; }
  virtual float get_best2_EMCal_4x4_eta() const { return 0; }
  virtual float get_best2_EMCal_4x4_phi() const { return 0; }

  // FullCalo 0.2x0.2
  virtual void set_best_FullCalo_0p2x0p2_E(const float /*E*/) {}
  virtual void set_best_FullCalo_0p2x0p2_eta(const float /*eta*/) {}
  virtual void set_best_FullCalo_0p2x0p2_phi(const float /*phi*/) {}

  virtual float get_best_FullCalo_0p2x0p2_E() const { return 0; }
  virtual float get_best_FullCalo_0p2x0p2_eta() const { return 0; }
  virtual float get_best_FullCalo_0p2x0p2_phi() const { return 0; }

  // FullCalo 0.4x0.4
  virtual void set_best_FullCalo_0p4x0p4_E(const float /*E*/) {}
  virtual void set_best_FullCalo_0p4x0p4_eta(const float /*eta*/) {}
  virtual void set_best_FullCalo_0p4x0p4_phi(const float /*phi*/) {}

  virtual float get_best_FullCalo_0p4x0p4_E() const { return 0; }
  virtual float get_best_FullCalo_0p4x0p4_eta() const { return 0; }
  virtual float get_best_FullCalo_0p4x0p4_phi() const { return 0; }

  // FullCalo 0.6x0.6
  virtual void set_best_FullCalo_0p6x0p6_E(const float /*E*/) {}
  virtual void set_best_FullCalo_0p6x0p6_eta(const float /*eta*/) {}
  virtual void set_best_FullCalo_0p6x0p6_phi(const float /*phi*/) {}

  virtual float get_best_FullCalo_0p6x0p6_E() const { return 0; }
  virtual float get_best_FullCalo_0p6x0p6_eta() const { return 0; }
  virtual float get_best_FullCalo_0p6x0p6_phi() const { return 0; }

  // FullCalo 0.8x0.8
  virtual void set_best_FullCalo_0p8x0p8_E(const float /*E*/) {}
  virtual void set_best_FullCalo_0p8x0p8_eta(const float /*eta*/) {}
  virtual void set_best_FullCalo_0p8x0p8_phi(const float /*phi*/) {}

  virtual float get_best_FullCalo_0p8x0p8_E() const { return 0; }
  virtual float get_best_FullCalo_0p8x0p8_eta() const { return 0; }
  virtual float get_best_FullCalo_0p8x0p8_phi() const { return 0; }

  // FullCalo 1.0x1.0
  virtual void set_best_FullCalo_1p0x1p0_E(const float /*E*/) {}
  virtual void set_best_FullCalo_1p0x1p0_eta(const float /*eta*/) {}
  virtual void set_best_FullCalo_1p0x1p0_phi(const float /*phi*/) {}

  virtual float get_best_FullCalo_1p0x1p0_E() const { return 0; }
  virtual float get_best_FullCalo_1p0x1p0_eta() const { return 0; }
  virtual float get_best_FullCalo_1p0x1p0_phi() const { return 0; }

 protected:
  CaloTriggerInfo() {}

 private:
  ClassDefOverride(CaloTriggerInfo, 1);
};

#endif  // TRIGGER_CALOTRIGGERINFO_H

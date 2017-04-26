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
  virtual void set_best_2x2_E(float E) {}
  virtual void set_best_2x2_eta(float eta) {}
  virtual void set_best_2x2_phi(float phi) {}
  virtual float get_best_2x2_E() { return 0; }
  virtual float get_best_2x2_eta() { return 0; }
  virtual float get_best_2x2_phi() { return 0; }
  virtual void set_best_4x4_E(float E) {}
  virtual void set_best_4x4_eta(float eta) {}
  virtual void set_best_4x4_phi(float phi) {}
  virtual float get_best_4x4_E() { return 0; }
  virtual float get_best_4x4_eta() { return 0; }
  virtual float get_best_4x4_phi() { return 0; }

 protected:
  CaloTriggerInfo() {}

 private:

  ClassDef(CaloTriggerInfo, 1);
};

#endif  // __CALOTRIGGERINFO_H__

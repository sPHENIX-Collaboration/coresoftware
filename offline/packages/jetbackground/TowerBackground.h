#ifndef JETBACKGROUND_TOWERBACKGROUND_H
#define JETBACKGROUND_TOWERBACKGROUND_H

#include <phool/PHObject.h>

#include <vector>

class TowerBackground : public PHObject
{
 public:
  ~TowerBackground() override {};

  void identify(std::ostream &os = std::cout) const override { os << "TowerBackground base class" << std::endl; };
  int isValid() const override { return 0; }

  virtual void set_UE(int /*layer*/, const std::vector<float> & /*UE*/) {}
  virtual void set_v2(float) {}
  virtual void set_Psi2(float) {}
  virtual void set_nStripsUsedForFlow(int) {}
  virtual void set_nTowersUsedForBkg(int) {}
  virtual void set_flow_failure_flag(bool) {}

  virtual std::vector<float> get_UE(int /*layer*/) const { return std::vector<float>(); };
  virtual float get_v2() const { return 0; }
  virtual float get_Psi2() const { return 0; }
  virtual int get_nStripsUsedForFlow() const { return 0; }
  virtual int get_nTowersUsedForBkg() const { return 0; }
  virtual bool get_flow_failure_flag() const { return false; }

 protected:
  TowerBackground() {}

 private:
  ClassDefOverride(TowerBackground, 1);
};

#endif

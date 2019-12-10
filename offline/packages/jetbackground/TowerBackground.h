#ifndef JETBACKGROUND_TOWERBACKGROUND_H
#define JETBACKGROUND_TOWERBACKGROUND_H

#include <phool/PHObject.h>

#include <vector>

class TowerBackground : public PHObject
{
 public:
  virtual ~TowerBackground(){};

  virtual void identify(std::ostream &os = std::cout) const { os << "TowerBackground base class" << std::endl; };
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  virtual void set_UE(int layer, const std::vector<float> &UE) {}
  virtual void set_v2(float v2) {}
  virtual void set_Psi2(float Psi2) {}
  virtual void set_nStripsUsedForFlow(int nStrips) {}
  virtual void set_nTowersUsedForBkg(int nTowers) {}

  virtual std::vector<float> get_UE(int layer) { return std::vector<float>(); };
  virtual float get_v2() { return 0; }
  virtual float get_Psi2() { return 0; }
  virtual int get_nStripsUsedForFlow() { return 0; }
  virtual int get_nTowersUsedForBkg() { return 0; }

 protected:
  TowerBackground() {}

 private:
  ClassDef(TowerBackground, 1);
};

#endif

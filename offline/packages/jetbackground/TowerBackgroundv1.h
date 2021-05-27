#ifndef JETBACKGROUND_TOWERBACKGROUNDV1_H
#define JETBACKGROUND_TOWERBACKGROUNDV1_H

#include "TowerBackground.h"

#include <iostream>
#include <vector>

class TowerBackgroundv1 : public TowerBackground
{
 public:
  TowerBackgroundv1();
  ~TowerBackgroundv1() override {}

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}
  int isValid() const override { return 1; }

  void set_UE(int layer, const std::vector<float> &UE) override { _UE[layer] = UE; }
  void set_v2(float v2) override { _v2 = v2; }
  void set_Psi2(float Psi2) override { _Psi2 = Psi2; }
  void set_nStripsUsedForFlow(int nStrips) override { _nStrips = nStrips; }
  void set_nTowersUsedForBkg(int nTowers) override { _nTowers = nTowers; }

  std::vector<float> get_UE(int layer) override { return _UE[layer]; }
  float get_v2() override { return _v2; }
  float get_Psi2() override { return _Psi2; }
  int get_nStripsUsedForFlow() override { return _nStrips; }

// our own - not from parent class
  virtual int get_nTowersUsedForFlow() { return _nTowers; }

 private:
  std::vector<std::vector<float> > _UE;
  float _v2;
  float _Psi2;
  int _nStrips;
  int _nTowers;

  ClassDefOverride(TowerBackgroundv1, 1);
};

#endif

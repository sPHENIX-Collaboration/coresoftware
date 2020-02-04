#ifndef JETBACKGROUND_TOWERBACKGROUNDV1_H
#define JETBACKGROUND_TOWERBACKGROUNDV1_H

#include "TowerBackground.h"

#include <iostream>
#include <vector>

class TowerBackgroundv1 : public TowerBackground
{
 public:
  TowerBackgroundv1();
  virtual ~TowerBackgroundv1() {}

  void identify(std::ostream &os = std::cout) const;
  void Reset() {}
  int isValid() const { return 1; }

  virtual void set_UE(int layer, const std::vector<float> &UE) { _UE[layer] = UE; }
  virtual void set_v2(float v2) { _v2 = v2; }
  virtual void set_Psi2(float Psi2) { _Psi2 = Psi2; }
  virtual void set_nStripsUsedForFlow(int nStrips) { _nStrips = nStrips; }
  virtual void set_nTowersUsedForBkg(int nTowers) { _nTowers = nTowers; }

  virtual std::vector<float> get_UE(int layer) { return _UE[layer]; }
  virtual float get_v2() { return _v2; }
  virtual float get_Psi2() { return _Psi2; }
  virtual int get_nStripsUsedForFlow() { return _nStrips; }
  virtual int get_nTowersUsedForFlow() { return _nTowers; }

 private:
  std::vector<std::vector<float> > _UE;
  float _v2;
  float _Psi2;
  int _nStrips;
  int _nTowers;

  ClassDef(TowerBackgroundv1, 1);
};

#endif

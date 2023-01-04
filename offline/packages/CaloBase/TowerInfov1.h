#ifndef TOWERINFOV1_H
#define TOWERINFOV1_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

#include <cmath>

class TowerInfov1 : public TowerInfo
{
 public:
  TowerInfov1() {}
  ~TowerInfov1() override {}
  void Reset() override;

  void setTime(short t) override { _time = t; }
  short getTime() override { return _time; }
  void setAmplitude(float amp) override { _amplitude = amp; }
  float getAmplitude() override { return _amplitude; }

 private:
  short _time = 0;
  float _amplitude = NAN;

  ClassDefOverride(TowerInfov1, 1);
};

#endif

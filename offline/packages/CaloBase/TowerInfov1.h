#ifndef TOWERINFOV1_H
#define TOWERINFOV1_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

class TowerInfov1 : public TowerInfo
{
 public:
  TowerInfov1();
  TowerInfov1(const TowerInfov1 &ti);
  ~TowerInfov1() override;
  void Reset() override;

  void setTime(short t) override { _time = t; }
  short getTime() override { return _time; }
  void setAmplitude(float amp) override { _amplitude = amp; }
  float getAmplitude() override { return _amplitude; }

private:
  short _time;
  float _amplitude;

  ClassDefOverride(TowerInfov1, 1);
};

#endif

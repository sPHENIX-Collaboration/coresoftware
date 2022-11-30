#ifndef TOWERINFO_H
#define TOWERINFO_H

#include <phool/PHObject.h>

#include <cmath>

class TowerInfo : public PHObject
{
 public:
  TowerInfo() = default;
  ~TowerInfo() override = default;
  void Reset() override;

  virtual void setTime(short /*t*/) { return; }
  virtual float getTime() { return NAN; }
  virtual void setAmplitude(float /*amp*/) { return; }
  virtual float getAmplitude() { return NAN; }

private:
  ClassDefOverride(TowerInfo, 1);
};

#endif

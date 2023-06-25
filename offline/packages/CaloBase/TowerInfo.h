#ifndef TOWERINFO_H
#define TOWERINFO_H

#include <phool/PHObject.h>

#include <cmath>

class TowerInfo : public PHObject
{
 public:
  TowerInfo() = default;
  ~TowerInfo() override = default;
  void Reset() override { return; }

  virtual void set_time(short /*t*/) { return; }
  virtual short get_time() { return -1; }
  virtual void set_energy(float /*energy*/) { return; }
  virtual float get_energy() { return NAN; }

 private:
  ClassDefOverride(TowerInfo, 1);
};

#endif

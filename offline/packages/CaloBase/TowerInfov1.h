#ifndef TOWERINFOV1_H
#define TOWERINFOV1_H

#include "TowerInfo.h"

class TowerInfov1 : public TowerInfo
{
 public:
  TowerInfov1() = default;
  TowerInfov1(TowerInfo& tower);
  ~TowerInfov1() override = default;
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  void set_time(float t) override { _time = t * 1000; }
  float get_time() override { return _time / 1000.; }
  void set_time_short(short t) override { _time = t * 1000; }
  short get_time_short() override { return short(_time / 1000); }
  void set_energy(float energy) override { _energy = energy; }
  float get_energy() override { return _energy; }
  void copy_tower(TowerInfo* tower) override;

 private:
  short _time{0};
  float _energy{0};

  ClassDefOverride(TowerInfov1, 1);
};

#endif

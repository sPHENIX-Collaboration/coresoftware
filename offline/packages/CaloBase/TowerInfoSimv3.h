#ifndef TOWERINFOSIMV3_H
#define TOWERINFOSIMV3_H

#include "TowerInfoSimv1.h"

#include <cstdint>  // For int16_t
#include <vector>

class TowerInfoSimv3 : public TowerInfoSimv1
{
 public:
  TowerInfoSimv3() = default;
  ~TowerInfoSimv3() override = default;

  void Reset() override;
  void Clear(Option_t* = "") override;

  void copy_tower(TowerInfo* tower) override;

  void set_nsample(int nsample) override;
  int get_nsample() const override { return _waveform.size(); }
  int16_t get_waveform_value(int index) const override;
  void set_waveform_value(int index, int16_t value) override;

 private:
  EdepMap _hitedeps;
  ShowerEdepMap _showeredeps;
  std::vector<int16_t> _waveform;

  ClassDefOverride(TowerInfoSimv3, 1);
  // Inherit other methods and properties from TowerInfoSimv1
};

#endif

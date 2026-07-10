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

  void copy_tower(TowerInfo* tower) override;

  void set_nsample(int nsample) override;
  int get_nsample() const override { return _waveform.size(); }
  int16_t get_waveform_value(int index) const override;
  void set_waveform_value(int index, int16_t value) override;

 private:
  std::vector<int16_t> _waveform;

  ClassDefOverride(TowerInfoSimv3, 1);
};

#endif

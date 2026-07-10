#ifndef TOWERINFOV5_H
#define TOWERINFOV5_H

#include "TowerInfov2.h"

#include <cstdint>  // For int16_t
#include <iosfwd>   // for ostream
#include <vector>

class TowerInfov5 : public TowerInfov2
{
 public:
  TowerInfov5() = default;
  ~TowerInfov5() override = default;

  void Reset() override;

  void identify(std::ostream& os) const override;

  void copy_tower(TowerInfo* tower) override;

  // Getter and setter for waveform
  void set_nsample(int nsample) override;
  int get_nsample() const override { return _waveform.size(); }
  int16_t get_waveform_value(int index) const override;
  void set_waveform_value(int index, int16_t value) override;

 private:
  std::vector<int16_t> _waveform;

  ClassDefOverride(TowerInfov5, 1);
  // Inherit other methods and properties from TowerInfov2
};

#endif

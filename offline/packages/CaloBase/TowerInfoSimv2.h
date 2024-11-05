#ifndef TOWERINFOSIMV2_H
#define TOWERINFOSIMV2_H

#include "TowerInfoSimv1.h"
#include "TowerInfo.h"

#include <cstdint>  // For int16_t

class TowerInfoSimv2 : public TowerInfoSimv1
{
 public:
    
  TowerInfoSimv2() {}
  ~TowerInfoSimv2() override {}
  
  void Reset() override;
  void Clear(Option_t* = "") override;

  void copy_tower(TowerInfo* tower) override;

  int get_nsample() const override { return nsample; }
  int16_t get_waveform_value(int index) const override;
  void set_waveform_value(int index, int16_t value) override;

 private:
  EdepMap _hitedeps;
  ShowerEdepMap _showeredeps;

  static const int nsample = 31;
  int16_t _waveform[nsample] = {0}; // Initializes the entire array to zero

  ClassDefOverride(TowerInfoSimv2, 1);
  // Inherit other methods and properties from TowerInfoSimv1
};

#endif

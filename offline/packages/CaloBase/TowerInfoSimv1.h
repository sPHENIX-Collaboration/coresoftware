#ifndef TOWERINFOSIMV1_H
#define TOWERINFOSIMV1_H

#include "TowerInfov2.h"
#include "TowerInfo.h"

#include <cstdint>  // For int16_t

class TowerInfoSimv1 : public TowerInfov2
{
 public:
    
  TowerInfoSimv1() {}
  ~TowerInfoSimv1() override {}
  
  void Reset() override;
  void Clear(Option_t* = "") override;

  void copy_tower(TowerInfo* tower) override;

  EdepMap& get_hitEdepMap() override;
  ShowerEdepMap& get_showerEdepMap() override;
  const EdepMap& get_hitEdepMap() const override;
  const ShowerEdepMap& get_showerEdepMap() const override;
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) override;
  void add_shower_edep(const int showerid, const float edep) override;

 private:
  EdepMap _hitedeps;
  ShowerEdepMap _showeredeps;

  ClassDefOverride(TowerInfoSimv1, 1);
  // Inherit other methods and properties from TowerInfov2
};

#endif

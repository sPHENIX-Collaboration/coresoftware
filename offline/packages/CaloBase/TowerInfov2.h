#ifndef TOWERINFOV2_H
#define TOWERINFOV2_H

#include "TowerInfo.h"
#include "TowerInfov1.h"

class TowerInfov2 : public TowerInfov1
{
 public:
  TowerInfov2() {}

  ~TowerInfov2() override {}

  void Reset() override;
  void Clear(Option_t* = "") override;
  void set_time(short t) override { TowerInfov1::set_time(t * 1000); }
  short get_time() override { return TowerInfov1::get_time() / 1000; }

  void set_time_float(float t) override { TowerInfov1::set_time(t * 1000); }
  float get_time_float() override { return TowerInfov1::get_time() / 1000.; }

  void set_chi2(float chi2) override { _chi2 = chi2; }
  float get_chi2() override { return _chi2; }
  void set_pedestal(float pedestal) override { _pedestal = pedestal; }
  float get_pedestal() override { return _pedestal; }

  void set_isHot(bool isHot) override { set_status_bit(0, isHot); }
  bool get_isHot() const override { return get_status_bit(0); }

  void set_isBadTime(bool isBadTime) override { set_status_bit(1, isBadTime); }
  bool get_isBadTime() const override { return get_status_bit(1); }

  void set_isBadChi2(bool isBadChi2) override { set_status_bit(2, isBadChi2); }
  bool get_isBadChi2() const override { return get_status_bit(2); }

  void set_isNotInstr(bool isNotInstr) override { set_status_bit(3, isNotInstr); }
  bool get_isNotInstr() const override { return get_status_bit(3); }

  void set_isNoCalib(bool isNoCalib) override { set_status_bit(4, isNoCalib); }
  bool get_isNoCalib() const override { return get_status_bit(4); }

  bool get_isGood() const override { return !((bool) _status); }

  uint8_t get_status() const override { return _status; }

  void set_status(uint8_t status) override { _status = status; }

  void copy_tower(TowerInfo* tower) override;

 private:
  float _chi2 = 0;
  float _pedestal = 0;
  uint8_t _status = 0;

  void set_status_bit(int bit, bool value)
  {
    if (bit < 0 || bit > 7)
    {
      return;
    }
    _status &= ~((uint8_t) 1 << bit);
    _status |= (uint8_t) value << bit;
  }

  bool get_status_bit(int bit) const
  {
    if (bit < 0 || bit > 7)
    {
      return false;  // default behavior
    }
    return (_status & ((uint8_t) 1 << bit)) != 0;
  }

  ClassDefOverride(TowerInfov2, 1);
  // Inherit other methods and properties from TowerInfov1
};

#endif

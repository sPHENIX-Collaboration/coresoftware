#ifndef TOWERINFOV2_H
#define TOWERINFOV2_H

#include "TowerInfov1.h"

class TowerInfov2 : public TowerInfov1
{
 public:
  TowerInfov2() {}
  TowerInfov2(TowerInfo& tower)
    : TowerInfov1(tower)
  {
  }
  ~TowerInfov2() override {}

  void Reset() override;
  void Clear(Option_t* = "") override;

  void set_time(short t) override { _time = t * 1000; }
  short get_time() override { return _time / 1000; }

  void set_time_float(float t) { _time = t * 1000; }
  float get_time_float() { return _time / 1000.; }

  void set_chi2(float chi2) { _chi2 = chi2; }
  float get_chi2() { return _chi2; }
  void set_pedestal(float pedestal) { _pedestal = pedestal; }
  float get_pedestal() { return _pedestal; }

  void set_isHot(bool isHot) { set_status_bit(0, isHot); }
  bool get_isHot() const { return get_status_bit(0); }

  void set_isBadTime(bool isBadTime) { set_status_bit(1, isBadTime); }
  bool get_isBadTime() const { return get_status_bit(1); }

  void set_isBadChi2(bool isBadChi2) { set_status_bit(2, isBadChi2); }
  bool get_isBadChi2() const { return get_status_bit(2); }

  void set_isNotInstr(bool isNotInstr) { set_status_bit(3, isNotInstr); }
  bool get_isNotInstr() const { return get_status_bit(3); }

  bool get_isGood() const { return !((bool) _status); }

  uint8_t get_status() const { return _status; }

  void set_status(uint8_t status) { _status = status; }

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

  // Inherit other methods and properties from TowerInfov1
};

#endif

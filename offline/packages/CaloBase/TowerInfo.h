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
  virtual void copy_tower(TowerInfo* /*tower*/) { return; }
  // methods in v2
  virtual void set_time_float(float /*t*/) { return; }
  virtual float get_time_float() { return NAN; }
  virtual void set_chi2(float /*chi2*/) { return; }
  virtual float get_chi2() { return NAN; }
  virtual void set_pedestal(float /*pedestal*/) { return; }
  virtual float get_pedestal() { return NAN; }
  virtual void set_isHot(bool /*isHot*/) { return; }
  virtual bool get_isHot() const { return false; }
  virtual void set_isBadTime(bool /*isBadTime*/) { return; }
  virtual bool get_isBadTime() const { return false; }
  virtual void set_isBadChi2(bool /*isBadChi2*/) { return; }
  virtual bool get_isBadChi2() const { return false; }
  virtual void set_isNotInstr(bool /*isNotInstr*/) { return; }
  virtual bool get_isNotInstr() const { return false; }
  virtual bool get_isGood() const { return true; }
  virtual uint8_t get_status() const { return 0; }
  virtual void set_status(uint8_t /*status*/) { return; }
  // methods in v3
  virtual int get_nsample() const {return 0;}
  virtual int16_t get_waveform_value(int /*index*/) const {return -1;}
  virtual void set_waveform_value(int /*index*/, int16_t /*value*/) {return;}

 private:
  ClassDefOverride(TowerInfo, 1);
};

#endif

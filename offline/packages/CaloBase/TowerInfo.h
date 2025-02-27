#ifndef TOWERINFO_H
#define TOWERINFO_H

#include <phool/PHObject.h>

#include <g4main/PHG4HitDefs.h>

#include <limits>
#include <map>

class TowerInfo : public PHObject
{
 public:
  typedef std::map<PHG4HitDefs::keytype, float> EdepMap;
  typedef std::map<int, float> ShowerEdepMap;
  TowerInfo() = default;
  ~TowerInfo() override = default;
  void Reset() override { return; }

  virtual void set_time(short /*t*/) { return; }
  virtual short get_time() { return -1; }
  virtual void set_energy(float /*energy*/) { return; }
  virtual float get_energy() { return std::numeric_limits<float>::signaling_NaN(); }
  virtual void copy_tower(TowerInfo* /*tower*/) { return; }
  // methods in v2
  virtual void set_time_float(float /*t*/) { return; }
  virtual float get_time_float() { return std::numeric_limits<float>::signaling_NaN(); }
  virtual void set_chi2(float /*chi2*/) { return; }
  virtual float get_chi2() { return std::numeric_limits<float>::signaling_NaN(); }
  virtual void set_pedestal(float /*pedestal*/) { return; }
  virtual float get_pedestal() { return std::numeric_limits<float>::signaling_NaN(); }
  virtual void set_isHot(bool /*isHot*/) { return; }
  virtual bool get_isHot() const { return false; }
  virtual void set_isBadTime(bool /*isBadTime*/) { return; }
  virtual bool get_isBadTime() const { return false; }
  virtual void set_isBadChi2(bool /*isBadChi2*/) { return; }
  virtual bool get_isBadChi2() const { return false; }
  virtual void set_isNotInstr(bool /*isNotInstr*/) { return; }
  virtual bool get_isNotInstr() const { return false; }
  virtual void set_isNoCalib(bool /*isNotInstr*/) { return; }
  virtual bool get_isNoCalib() const { return false; }
  virtual void set_isZS(bool /*isZS*/) { return; }
  virtual bool get_isZS() const { return false; }
  virtual void set_isRecovered(bool /*isRecovered*/) { return; }
  virtual bool get_isRecovered() const { return false; }
  virtual void set_isSaturated(bool /*isSaturated*/) { return; }
  virtual bool get_isSaturated() const { return false; }
  virtual bool get_isGood() const { return true; }
  virtual uint8_t get_status() const { return 0; }
  virtual void set_status(uint8_t /*status*/) { return; }
  // methods in v3
  virtual int get_nsample() const { return 0; }
  virtual int16_t get_waveform_value(int /*index*/) const { return -1; }
  virtual void set_waveform_value(int /*index*/, int16_t /*value*/) { return; }
  // methods in sim v1
  virtual EdepMap& get_hitEdepMap() { static EdepMap dummy; return dummy; }
  virtual ShowerEdepMap& get_showerEdepMap() { static ShowerEdepMap dummy; return dummy; }
  virtual const EdepMap& get_hitEdepMap() const { static EdepMap dummy; return dummy; }
  virtual const ShowerEdepMap& get_showerEdepMap() const { static ShowerEdepMap dummy; return dummy; }
  virtual void add_edep(const PHG4HitDefs::keytype /*g4hitid*/, const float /*edep*/) { return; }
  virtual void add_shower_edep(const int /*showerid*/, const float /*edep*/) { return; }

 private:
  ClassDefOverride(TowerInfo, 1);
};

#endif

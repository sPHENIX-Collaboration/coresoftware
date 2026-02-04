#ifndef SEPD_EVENTPLANECALIB_EVENTPLANEDATA_H
#define SEPD_EVENTPLANECALIB_EVENTPLANEDATA_H

#include <phool/PHObject.h>

#include <array>
#include <limits>

class EventPlaneData : public PHObject
{
 public:
  EventPlaneData();
  ~EventPlaneData() override = default;

  void Reset() override {*this = EventPlaneData();} // check if this works
  // this should be in an sepd define (e.g. ../../../offline/packages/epd/EPDDefs.h)
  static constexpr int SEPD_CHANNELS = 744;
  void set_event_id(int id) {event_id = id;}
  int get_event_id() const {return event_id;}

  void set_event_zvertex(double vtx) {event_zvertex = vtx;}
  double get_event_zvertex() const {return event_zvertex;}

  void set_sepd_totalcharge(double chg) {sepd_totalcharge = chg;}
  double get_sepd_totalcharge() const {return sepd_totalcharge;}

  void set_sepd_charge(int channel, double chg) {sepd_charge[channel] = chg;}
  double get_sepd_charge(int channel) const {return sepd_charge[channel];}

  void set_sepd_phi(int channel, double phi) {sepd_phi[channel] = phi;}
  double get_sepd_phi(int channel) const {return sepd_phi[channel];}
  
 private:
  int event_id {0};
  double event_zvertex {std::numeric_limits<double>::quiet_NaN()};
  double event_centrality{std::numeric_limits<double>::quiet_NaN()};
  double sepd_totalcharge{std::numeric_limits<double>::quiet_NaN()};
  
  std::array<double, SEPD_CHANNELS> sepd_charge {};
  std::array<double, SEPD_CHANNELS> sepd_phi {};
  ClassDefOverride(EventPlaneData, 1);
};

#endif

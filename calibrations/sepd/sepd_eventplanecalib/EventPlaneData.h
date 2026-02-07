#ifndef SEPD_EVENTPLANECALIB_EVENTPLANEDATA_H
#define SEPD_EVENTPLANECALIB_EVENTPLANEDATA_H

#include "QVecDefs.h"

#include <phool/PHObject.h>

#include <array>
#include <cmath>
#include <limits>

class EventPlaneData : public PHObject
{
 public:
  EventPlaneData();
  ~EventPlaneData() override = default;

  EventPlaneData(const EventPlaneData&) = default;
  EventPlaneData& operator=(const EventPlaneData&) = default;
  EventPlaneData(EventPlaneData&&) = default;
  EventPlaneData& operator=(EventPlaneData&&) = default;

  void Reset() override;
  void set_event_id(int id) {event_id = id;}
  int get_event_id() const {return event_id;}

  void set_event_zvertex(double vtx) {event_zvertex = vtx;}
  double get_event_zvertex() const {return event_zvertex;}

  void set_sepd_totalcharge(double chg) {sepd_totalcharge = chg;}
  double get_sepd_totalcharge() const {return sepd_totalcharge;}

  void set_sepd_charge(int channel, double chg) {sepd_charge[channel] = chg;}
  double get_sepd_charge(int channel) const {return sepd_charge[channel];}

  void set_event_centrality(double cent) { event_centrality = cent; }
  double get_event_centrality() const { return event_centrality; }

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;

 private:
  int event_id {0};
  double event_zvertex {std::numeric_limits<double>::quiet_NaN()};
  double event_centrality{std::numeric_limits<double>::quiet_NaN()};
  double sepd_totalcharge{std::numeric_limits<double>::quiet_NaN()};
  
  std::array<double, QVecShared::SEPD_CHANNELS> sepd_charge {};
  ClassDefOverride(EventPlaneData, 1);
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICEVENTHEADERV1_H
#define EICEVENTHEADERV1_H

#include "EicEventHeader.h"

#include <cstdint>
#include <iostream>
#include <map>

class EicEventHeaderv1 : public EicEventHeader
{
 public:
  EicEventHeaderv1() {}
  explicit EicEventHeaderv1(const EicEventHeader *eicevt);
  ~EicEventHeaderv1() override {}

  //  void identify(std::ostream& os  = std::cout) const;
  void Reset() override;

  // property relates methods
  bool has_property(const PROPERTY prop_id) const override;
  float get_property_float(const PROPERTY prop_id) const override;
  int get_property_int(const PROPERTY prop_id) const override;
  unsigned int get_property_uint(const PROPERTY prop_id) const override;
  void set_property(const PROPERTY prop_id, const float value) override;
  void set_property(const PROPERTY prop_id, const int value) override;
  void set_property(const PROPERTY prop_id, const unsigned int value) override;

  // Generator specific values
  void set_eventgenerator_type(const int i) override { set_property(prop_eventgen, i); }
  int get_eventgenerator_type() const override { return get_property_int(prop_eventgen); }

  // Milou
  void set_milou_weight(const float val) override { set_property(prop_milou_weight, val); }
  float get_milou_weight() const override { return get_property_float(prop_milou_weight); }
  void set_milou_trueX(const float val) override { set_property(prop_milou_truex, val); }
  float get_milou_trueX() const override { return get_property_float(prop_milou_truex); }
  void set_milou_trueQ2(const float val) override { set_property(prop_milou_trueq2, val); }
  float get_milou_trueQ2() const override { return get_property_float(prop_milou_trueq2); }

  // DEMP
  void set_demp_weight(const float val) override { set_property(prop_demp_weight, val); }
  float get_demp_weight() const override { return get_property_float(prop_demp_weight); }

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const override;
  void set_property_nocheck(const PROPERTY prop_id, const unsigned int ui) override { prop_map[prop_id] = ui; }

  //! storage types for properties
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;
  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property {
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in)
      : idata(in)
    {
    }
    u_property(uint32_t in)
      : uidata(in)
    {
    }
    u_property(float in)
      : fdata(in)
    {
    }
    u_property()
      : uidata(0)
    {
    }

    prop_storage_t get_storage() const { return uidata; }
  };

  //! container for  properties
  prop_map_t prop_map;

  ClassDefOverride(EicEventHeaderv1, 1)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICEVENTHEADER_H
#define EICEVENTHEADER_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

class EicEventHeader : public PHObject
{
 public:
  EicEventHeader() {}
  ~EicEventHeader() override;

  void identify(std::ostream& os = std::cout) const override;
  void CopyFrom(const PHObject* phobj) override;

  void Reset() override;

  virtual void set_eventgenerator_type(const int) { return; }
  virtual int get_eventgenerator_type() const { return -99999; }

  // Milou

  virtual void set_milou_weight(const float) { return; }
  virtual float get_milou_weight() const { return NAN; }
  virtual void set_milou_trueX(const float) { return; }
  virtual float get_milou_trueX() const { return NAN; }
  virtual void set_milou_trueQ2(const float) { return; }
  virtual float get_milou_trueQ2() const { return NAN; }


//  void set_milou_weight(const float val) override { set_property(prop_milou_weight, val); }
//  float get_milou_weight() const override { return get_property_float(prop_milou_weight); }

  // DEMP
  virtual void set_demp_weight(const float) { return; }
  virtual float get_demp_weight() const { return NAN; }

  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to EicEventHeader::get_property_info
  enum PROPERTY
  {  //
    prop_eventgen = 1,
    prop_milou_weight = 2,
    prop_milou_truex = 3,
    prop_milou_trueq2 = 4,

    prop_demp_weight = 5,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE
  {  //
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  enum EvtGen
  {
    Milou = 1,
    DEMP = 2
  };

  virtual bool has_property(const PROPERTY /*prop_id*/) const { return false; }
  virtual float get_property_float(const PROPERTY /*prop_id*/) const { return NAN; }
  virtual int get_property_int(const PROPERTY /*prop_id*/) const { return INT_MIN; }
  virtual unsigned int get_property_uint(const PROPERTY /*prop_id*/) const { return UINT_MAX; }
  virtual void set_property(const PROPERTY /*prop_id*/, const float /*value*/) { return; }
  virtual void set_property(const PROPERTY /*prop_id*/, const int /*value*/) { return; }
  virtual void set_property(const PROPERTY /*prop_id*/, const unsigned int /*value*/) { return; }
  static std::pair<const std::string, PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:
  virtual unsigned int get_property_nocheck(const PROPERTY /*prop_id*/) const { return UINT_MAX; }
  virtual void set_property_nocheck(const PROPERTY /*prop_id*/, const unsigned int) { return; }

  std::map<std::string, double> evInfo;

  ClassDefOverride(EicEventHeader, 1)
};

#endif

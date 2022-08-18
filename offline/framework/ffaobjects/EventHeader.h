// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_EVENTHEADER_H
#define FFAOBJECTS_EVENTHEADER_H

#include <phool/PHObject.h>

#include <cmath>
#include <cstdint>  // for int64_t
#include <ctime>
#include <iostream>  // for cout, ostream
#include <string>    // for string

//! base class for EventHeaders
class EventHeader : public PHObject
{
 public:
  /// dtor
  ~EventHeader() override = default;

  /// Clear Event
  void Reset() override;

  /*
   * identify Function from PHObject
   * @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  /// get Run Number
  virtual int get_RunNumber() const { return 0; }
  /// set Run Number
  virtual void set_RunNumber(const int /*run*/) { return; }

  /// get Event Number
  virtual int get_EvtSequence() const { return 0; }
  /// set Event Number
  virtual void set_EvtSequence(const int /*ival*/) { return; }

  //! bunch crossing
  virtual void set_BunchCrossing(int64_t bcr) { set_intval("bcr", bcr); }

  //! bunch crossing
  virtual int64_t get_BunchCrossing() const { return get_intval("bcr"); }

  virtual void set_floatval(const std::string & /*name*/, const float /*fval*/) { return; }
  virtual float get_floatval(const std::string & /*name*/) const { return NAN; }

  virtual void set_intval(const std::string & /*name*/, const int64_t /*ival*/) { return; }
  virtual int64_t get_intval(const std::string & /*name*/) const { return -999999; }

  /// get Event Type (Data,rejected,EOR,BOR,...)
  int get_EvtType() const { return get_intval("type"); }
  /// set Event Type (Data,rejected,EOR,BOR,...)
  void set_EvtType(const int ival) { set_intval("type", ival); }

  void set_ImpactParameter(const double rval) { set_floatval("bimp", rval); }
  float get_ImpactParameter() const { return get_floatval("bimp"); }

  void set_EventPlaneAngle(const double rval) { set_floatval("rplane", rval); }
  float get_EventPlaneAngle() const { return get_floatval("rplane"); }

  void set_eccentricity(const double rval) { set_floatval("ecc", rval); }
  float get_eccentricity() const { return get_floatval("ecc"); }

  void set_ncoll(const int ival) { set_intval("ncoll", ival); }
  int get_ncoll() const { return get_intval("ncoll"); }

  void set_npart(const int ival) { set_intval("npart", ival); }
  int get_npart() const { return get_intval("npart"); }

  void set_TimeStamp(const time_t tval) { set_intval("time", tval); }
  time_t get_TimeStamp() const { return get_intval("time"); }

 private:  // prevent doc++ from showing ClassDefOverride
  ClassDefOverride(EventHeader, 1)
};

#endif

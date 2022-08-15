// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_EVENTHEADERV1_H
#define FFAOBJECTS_EVENTHEADERV1_H

#include "EventHeader.h"

#include <cstdint>   // for int64_t
#include <iostream>  // for cout, ostream
#include <map>
#include <string>

class PHObject;

//! simple event header with ID and time
class EventHeaderv1 : public EventHeader
{
 public:
  /// ctor
  EventHeaderv1() = default;
  /// dtor
  ~EventHeaderv1() override = default;

  PHObject *CloneMe() const override { return new EventHeaderv1(*this); }

  ///  Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  /// get Run Number
  int get_RunNumber() const override { return RunNumber; }
  /// set Run Number
  void set_RunNumber(const int run) override { RunNumber = run; }

  /// get Event Number
  int get_EvtSequence() const override { return EvtSequence; }
  /// set Event Number
  void set_EvtSequence(const int evtno) override { EvtSequence = evtno; }

  void set_floatval(const std::string &name, const float fval) override;
  float get_floatval(const std::string &name) const override;

  void set_intval(const std::string &name, const int64_t ival) override;
  int64_t get_intval(const std::string &name) const override;

 private:
  int RunNumber = 0;    // Run number
  int EvtSequence = 0;  // Event number
  std::map<std::string, int64_t> m_IntEventProperties;
  std::map<std::string, float> m_FloatEventProperties;

  ClassDefOverride(EventHeaderv1, 2)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTHEADERV1_H
#define G4MAIN_PHG4EVENTHEADERV1_H

#include "PHG4EventHeader.h"

#include <iostream>

///
class PHG4EventHeaderv1: public PHG4EventHeader
{
 public:

  PHG4EventHeaderv1();

  /// dtor
  ~PHG4EventHeaderv1() override {}

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  /// get Event Number
  int get_EvtSequence() const override {return evtseq;}
  /// set Event Number
  void set_EvtSequence(const int ival) override {evtseq = ival;}

  float get_ImpactParameter() const override {return bimp;}
  void set_ImpactParameter(const float b) override {bimp = b;}

  float get_EventPlaneAngle() const override {return rplane;}
  void set_EventPlaneAngle(const float r) override {rplane = r;}

 protected:
  int evtseq;
  float bimp;
  float rplane;

 private: // prevent doc++ from showing ClassDefOverride
  ClassDefOverride(PHG4EventHeaderv1,1)

};

#endif

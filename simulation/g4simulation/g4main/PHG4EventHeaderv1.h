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
  virtual ~PHG4EventHeaderv1() {}

  /// Clear Event
  virtual void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const;

  /// isValid returns non zero if object contains valid data
  int isValid() const;

  /// get Event Number
  int get_EvtSequence() const {return evtseq;}
  /// set Event Number
  void set_EvtSequence(const int ival) {evtseq = ival;}

  float get_ImpactParameter() const {return bimp;}
  void set_ImpactParameter(const float b) {bimp = b;}

  float get_EventPlaneAngle() const {return rplane;}
  void set_EventPlaneAngle(const float r) {rplane = r;}

 protected:
  int evtseq;
  float bimp;
  float rplane;

 private: // prevent doc++ from showing ClassDef
  ClassDef(PHG4EventHeaderv1,1)

};

#endif

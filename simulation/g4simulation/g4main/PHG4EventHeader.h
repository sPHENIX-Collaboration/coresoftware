// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTHEADER_H
#define G4MAIN_PHG4EVENTHEADER_H

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>

///
class PHG4EventHeader: public PHObject
{
 public:

  /// dtor
  ~PHG4EventHeader() override {}

  /// Clear Event
  void Reset() override
    {
      std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
      return;
    }

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const override
    {
      os << "identify yourself: virtual PHG4EventHeader Object" << std::endl;
      return;
    }

  /// isValid returns non zero if object contains valid data
  int isValid() const override
    {
      std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
      return 0;
    }

  /// get Event Number
  virtual int get_EvtSequence() const {return -9999;}
  /// set Event Number
  virtual void set_EvtSequence(const int /*ival*/) {return;}

  virtual float get_ImpactParameter() const {return NAN;}
  virtual void set_ImpactParameter(const float) {return;}

  virtual float get_EventPlaneAngle() const {return NAN;}
  virtual void set_EventPlaneAngle(const float) {return;}


 private: // prevent doc++ from showing ClassDefOverride
  ClassDefOverride(PHG4EventHeader,1)

};

#endif




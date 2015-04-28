#ifndef __PHG4Hitv4_H__
#define __PHG4Hitv4_H__

#include "PHG4Hitv1.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv4;
extern G4Allocator<PHG4Hitv4> PHG4Hitv4Allocator;
#endif

class PHG4Hitv4 : public PHG4Hitv1
{
 public:
  PHG4Hitv4();
  PHG4Hitv4(PHG4Hit const &g4hit);
  // The indices here represent the entry and exit points of the particle
  float get_eion() const {return eion;}
  void set_eion(const float f) {eion = f;}

  virtual void print() const;

 protected:
  float eion;

  ClassDef(PHG4Hitv4,1)
};

#endif

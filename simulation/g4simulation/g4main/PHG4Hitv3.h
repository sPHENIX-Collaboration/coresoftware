#ifndef __PHG4Hitv3_H__
#define __PHG4Hitv3_H__

#include "PHG4Hitv1.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv3;
extern G4Allocator<PHG4Hitv3> PHG4Hitv3Allocator;
#endif

class PHG4Hitv3 : public PHG4Hitv1
{
 public:
  PHG4Hitv3();
  PHG4Hitv3(PHG4Hit const &g4hit);

  int get_scint_id() const {return scint_id;}

  void set_scint_id(const int i) {scint_id = i;}

  virtual void print() const;  

 protected:

  int scint_id;

  ClassDef(PHG4Hitv3,1)
};

#endif

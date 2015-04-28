#ifndef PHG4Hitv5_H__
#define PHG4Hitv5_H__

#include "PHG4Hitv4.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv5;
extern G4Allocator<PHG4Hitv5> PHG4Hitv5Allocator;
#endif

class PHG4Hitv5 : public PHG4Hitv4
{
 public:
  PHG4Hitv5();
  PHG4Hitv5(PHG4Hit const &g4hit);

  int get_scint_id() const {return scint_id;}

  void set_scint_id(const int i) {scint_id = i;}

  virtual void print() const;  

 protected:

  int scint_id;

  ClassDef(PHG4Hitv5,1)
};

#endif

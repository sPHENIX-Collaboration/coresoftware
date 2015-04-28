#ifndef __PHG4Hitv2_H__
#define __PHG4Hitv2_H__

#include "PHG4Hitv1.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv2;
extern G4Allocator<PHG4Hitv2> PHG4Hitv2Allocator;
#endif

class PHG4Hitv2 : public PHG4Hitv1
{
 public:
  PHG4Hitv2();
  PHG4Hitv2(PHG4Hit const &g4hit);
  // The indices here represent the entry and exit points of the particle
  float get_px(const int i) const {return px[i];}
  float get_py(const int i) const {return py[i];}
  float get_pz(const int i) const {return pz[i];}
  
  void set_px(const int i, const float f) {px[i]=f;}
  void set_py(const int i, const float f) {py[i]=f;}
  void set_pz(const int i, const float f) {pz[i]=f;}

  virtual void print() const;  

 protected:
  // Store both the entry and exit momentum of the particle
  // Remember, particles do not always enter on the inner edge!
  float px[2];
  float py[2];
  float pz[2];

  ClassDef(PHG4Hitv2,1)
};

#endif

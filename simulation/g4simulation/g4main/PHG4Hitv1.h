#ifndef __PHG4Hitv1_H__
#define __PHG4Hitv1_H__

#include "PHG4Hit.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv1;
extern G4Allocator<PHG4Hitv1> PHG4Hitv1Allocator;
#endif

class PHG4Hitv1 : public PHG4Hit
{
 public:
  PHG4Hitv1();
  PHG4Hitv1(PHG4Hit const &g4hit);
  // The indices here represent the entry and exit points of the particle
  float get_x(const int i) const {return x[i];}
  float get_y(const int i) const {return y[i];}
  float get_z(const int i) const {return z[i];}
  float get_t(const int i) const {return t[i];}
  float get_edep() const {return edep;}
  unsigned int get_layer() const {return layer;}
  unsigned int get_hit_id() const {return hitid;}
  int get_trkid() const {return trackid;}
  
  void set_x(const int i, const float f) {x[i]=f;}
  void set_y(const int i, const float f) {y[i]=f;}
  void set_z(const int i, const float f) {z[i]=f;}
  void set_t(const int i, const float f) {t[i]=f;}
  void set_edep(const float f) {edep = f;}
  void set_layer(const unsigned int i) {layer=i;}
  void set_hit_id(const unsigned int i) {hitid=i;}
  void set_trkid(const int i) {trackid=i;}

  virtual void print() const;

 protected:
  // Store both the entry and exit points of the particle
  // Remember, particles do not always enter on the inner edge!
  float x[2];
  float y[2];
  float z[2];
  float t[2];
  unsigned int layer;
  unsigned int hitid;
  int trackid;
  float edep;

  ClassDef(PHG4Hitv1,1)
};

#endif

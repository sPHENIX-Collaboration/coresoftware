#ifndef __PHG4Hitv7_H__
#define __PHG4Hitv7_H__

#include "PHG4Hitv2.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv7;
extern G4Allocator<PHG4Hitv7> PHG4Hitv7Allocator;
#endif

class PHG4Hitv7 : public PHG4Hitv2
{
 public:
  PHG4Hitv7();
  PHG4Hitv7(PHG4Hit const &g4hit);
  
  int get_strip_z_index() const {return strip_z_index;}
  int get_strip_y_index() const {return strip_y_index;}
  int get_ladder_z_index() const {return ladder_z_index;}
  int get_ladder_phi_index() const {return ladder_phi_index;}
  
  void set_strip_z_index(const int i) {strip_z_index = i;}
  void set_strip_y_index(const int i) {strip_y_index = i;}
  void set_ladder_z_index(const int i) {ladder_z_index = i;}
  void set_ladder_phi_index(const int i) {ladder_phi_index = i;}
  
  virtual void print() const;  

 protected:
  // Store the index numbers neede to find the location of the sensor strip that was hit
  // These give the location of the strip within its sensor
  int strip_z_index;
  int strip_y_index;
  // these give the location of the sensor within the tracker
  int ladder_z_index;
  int ladder_phi_index;
  
  ClassDef(PHG4Hitv7,1)
    };

#endif

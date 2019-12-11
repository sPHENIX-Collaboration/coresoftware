// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV2_H
#define G4DETECTORS_PHG4CYLINDERCELLV2_H

#include "PHG4CylinderCellv1.h"

#include <iostream>
#include <string>                // for string

class PHG4CylinderCellv2 : public PHG4CylinderCellv1
{
 public:

  PHG4CylinderCellv2();
  virtual ~PHG4CylinderCellv2(){}

  void identify(std::ostream& os = std::cout) const;

  void set_sensor_index(const std::string &si) {sensor_index = si;}
  std::string get_sensor_index() const  {return sensor_index;}

  void set_ladder_phi_index(const int i) {ladder_phi_index = i;}
  int get_ladder_phi_index() const {return ladder_phi_index;}

  void set_ladder_z_index(const int i) {ladder_z_index = i;}
  int get_ladder_z_index() const {return ladder_z_index;}
  
 protected:

  int ladder_phi_index;
  int ladder_z_index; 
  std::string sensor_index;
  
  ClassDef(PHG4CylinderCellv2,1)
};

#endif

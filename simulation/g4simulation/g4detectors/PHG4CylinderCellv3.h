#ifndef G4DETECTORS_PHG4CYLINDERCELLV3_H
#define G4DETECTORS_PHG4CYLINDERCELLV3_H

#include "PHG4CylinderCellv1.h"

#include <iostream>

class PHG4CylinderCellv3 : public PHG4CylinderCellv1
{
 public:

  PHG4CylinderCellv3();
  virtual ~PHG4CylinderCellv3(){}

  void identify(std::ostream& os = std::cout) const;

  void set_j_index(const int i) {j_index = i;}
  int get_j_index() const {return j_index;}

  void set_k_index(const int i) {k_index = i;}
  int get_k_index() const {return k_index;}

  void set_l_index(const int i) {l_index = i;}
  int get_l_index() const {return l_index;}
  
 protected:

  int j_index;
  int k_index; 
  int l_index; 
  
  ClassDef(PHG4CylinderCellv3,1)
};

#endif

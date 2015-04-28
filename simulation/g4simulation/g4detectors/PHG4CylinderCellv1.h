#ifndef PHG4CYLINDERCELLV1_H
#define PHG4CYLINDERCELLV1_H

#include "PHG4CylinderCell.h"
#include <cmath>
#include <map>

class PHG4CylinderCellv1 : public PHG4CylinderCell
{
 public:

  PHG4CylinderCellv1();
  virtual ~PHG4CylinderCellv1(){}

  void identify(std::ostream& os = std::cout) const;

  std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator > get_g4hits()
    {return make_pair(edeps.begin(), edeps.end());}
  void add_edep(const unsigned int g4hitid, const float edep);
  void set_cell_id(const unsigned int id) {cellid = id;}
  void set_layer(const unsigned int i) {layer = i;}
  double get_edep() const;
  unsigned int get_layer() const {return layer;}
  unsigned int get_cell_id() const {return cellid;}
  int get_binz() const {return binz;}
  int get_binphi() const {return binphi;}
  int get_bineta() const {return get_binz();}

  void set_zbin(const int i) {binz = i;}
  void set_etabin(const int i) {set_zbin(i);}
  void set_phibin(const int i) {binphi = i;}

 protected:

  unsigned int layer;
  unsigned int cellid;
  int binz;
  int binphi;
  std::map<unsigned int, float> edeps;
   
  ClassDef(PHG4CylinderCellv1,1)
};

#endif

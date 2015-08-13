#ifndef PHG4CYLINDERCELLV1_H
#define PHG4CYLINDERCELLV1_H

#include "PHG4CylinderCell.h"
#include <g4main/PHG4HitDefs.h>

#include <cmath>
#include <map>

class PHG4CylinderCellv1 : public PHG4CylinderCell
{
 public:

  PHG4CylinderCellv1();
  virtual ~PHG4CylinderCellv1(){}

  void identify(std::ostream& os = std::cout) const;

   EdepConstRange get_g4hits()
    {return make_pair(edeps.begin(), edeps.end());}
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);
  void set_cell_id(const PHG4CylinderCellDefs::keytype id) {cellid = id;}
  void set_layer(const unsigned int i) {layer = i;}
  double get_edep() const;
  unsigned int get_layer() const {return layer;}
  PHG4CylinderCellDefs::keytype get_cell_id() const {return cellid;}
  int get_binz() const {return binz;}
  int get_binphi() const {return binphi;}
  int get_bineta() const {return get_binz();}

  void set_zbin(const int i) {binz = i;}
  void set_etabin(const int i) {set_zbin(i);}
  void set_phibin(const int i) {binphi = i;}

 protected:

  unsigned int layer;
  PHG4CylinderCellDefs::keytype cellid;
  int binz;
  int binphi;
  EdepMap edeps;
   
  ClassDef(PHG4CylinderCellv1,1)
};

#endif

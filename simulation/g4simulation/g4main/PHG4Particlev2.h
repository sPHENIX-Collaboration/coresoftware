#ifndef PHG4PARTICLEV2_H__
#define PHG4PARTICLEV2_H__

#include "PHG4Particlev1.h"

class PHG4Particlev2: public PHG4Particlev1
{
 public:
  PHG4Particlev2();
  PHG4Particlev2(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev2(const PHG4Particle *in);

  virtual ~PHG4Particlev2() {}

  int get_track_id() const {return trkid;}
  int get_vtx_id() const {return vtxid;}
  int get_parent_id() const {return parentid;}
  int get_primary_id() const {return primaryid;}
  double get_e() const {return fe;}

  // this function was implimented in PHG4Particlev1
//  int get_barcode() const {return barcode;}

  void set_track_id(const int i) {trkid = i;}
  void set_vtx_id(const int i) {vtxid = i;}
  void set_parent_id(const int i) {parentid = i;}
  void set_primary_id(const int i) {primaryid = i;}
  void set_e(const double e) {fe = e;}

  // this function was implimented in PHG4Particlev1
//  void set_barcode(const int bcd) {barcode = bcd;}

  void identify(std::ostream& os = std::cout) const;

 protected:
  int trkid;
  int vtxid;
  int parentid;
  int primaryid;
  double fe;
  // this variable was implimented in PHG4Particlev1
//  int barcode;

  ClassDef(PHG4Particlev2,2)
};


#endif

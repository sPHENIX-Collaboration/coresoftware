#ifndef G4MAIN_PHG4VTXPOINTV1_H
#define G4MAIN_PHG4VTXPOINTV1_H

#include "PHG4VtxPoint.h"

#include <climits> // for INT_MIN
#include <cmath> // def of NAN
#include <iostream>   // for cout, ostream

class PHG4VtxPointv1: public PHG4VtxPoint
{
 public:
  PHG4VtxPointv1():
    vx(NAN),
    vy(NAN),
    vz(NAN),
    t0(NAN),
    id(INT_MIN)
      {}

  PHG4VtxPointv1(const PHG4VtxPoint *vtx):
    vx(vtx->get_x()),
    vy(vtx->get_y()),
    vz(vtx->get_z()),
    t0(vtx->get_t()),
    id(vtx->get_id())
      {}

  PHG4VtxPointv1(const double x, const double y, const double z, const double t, const int id_value = INT_MIN):
    vx(x),
    vy(y),
    vz(z),
    t0(t),
    id(id_value)
      {}

  virtual ~PHG4VtxPointv1() {}

  void set_x(const double r) {vx = r;}
  void set_y(const double r) {vy = r;}
  void set_z(const double r) {vz = r;}
  void set_t(const double r) {t0 = r;}
  void set_id(const int i) {id = i;}

  double get_x() const {return vx;}
  double get_y() const {return vy;}
  double get_z() const {return vz;}
  double get_t() const {return t0;}
  int get_id() const {return id;}

  void identify(std::ostream& os = std::cout) const;

 protected:

  double vx;
  double vy;
  double vz;
  double t0;

  //! id tag for this vertex
  int id;

  ClassDef(PHG4VtxPointv1,2)
};




#endif

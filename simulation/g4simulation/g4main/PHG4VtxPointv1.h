#ifndef PHG4VTXPOINTV1_H__
#define PHG4VTXPOINTV1_H__

#include "PHG4VtxPoint.h"
#include <cmath> // def of NAN

class PHG4VtxPointv1: public PHG4VtxPoint
{
 public:
  PHG4VtxPointv1():
    vx(NAN),
    vy(NAN),
    vz(NAN),
    t0(NAN)
      {}

  PHG4VtxPointv1(const PHG4VtxPoint *vtx):
    vx(vtx->get_x()),
    vy(vtx->get_y()),
    vz(vtx->get_z()),
    t0(vtx->get_t())
      {}

  PHG4VtxPointv1(const double x, const double y, const double z, const double t):
    vx(x),
    vy(y),
    vz(z),
    t0(t)
      {}

  virtual ~PHG4VtxPointv1() {}

  void set_x(const double r) {vx = r;}
  void set_y(const double r) {vy = r;}
  void set_z(const double r) {vz = r;}
  void set_t(const double r) {t0 = r;}

  double get_x() const {return vx;}
  double get_y() const {return vy;}
  double get_z() const {return vz;}
  double get_t() const {return t0;}

  void identify(std::ostream& os = std::cout) const;

 protected:
  double vx;
  double vy;
  double vz;
  double t0;
  ClassDef(PHG4VtxPointv1,1)
};




#endif

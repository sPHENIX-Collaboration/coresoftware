// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINTV1_H
#define G4MAIN_PHG4VTXPOINTV1_H

#include "PHG4VtxPoint.h"

#include <climits>   // for INT_MIN
#include <cmath>     // def of NAN
#include <iostream>  // for cout, ostream

class PHG4VtxPointv1 : public PHG4VtxPoint
{
 public:
  PHG4VtxPointv1() = default;

  PHG4VtxPointv1(const PHG4VtxPoint* vtx)
    : vx(vtx->get_x())
    , vy(vtx->get_y())
    , vz(vtx->get_z())
    , t0(vtx->get_t())
    , id(vtx->get_id())
  {
  }

  PHG4VtxPointv1(const double x, const double y, const double z, const double t, const int id_value = std::numeric_limits<int>::min())
    : vx(x)
    , vy(y)
    , vz(z)
    , t0(t)
    , id(id_value)
  {
  }

  ~PHG4VtxPointv1() override = default;

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_x(const double r) override { vx = r; }
  void set_y(const double r) override { vy = r; }
  void set_z(const double r) override { vz = r; }
  void set_t(const double r) override { t0 = r; }
  void set_id(const int i) override { id = i; }

  double get_x() const override { return vx; }
  double get_y() const override { return vy; }
  double get_z() const override { return vz; }
  double get_t() const override { return t0; }
  int get_id() const override { return id; }

 protected:
  double vx{std::numeric_limits<double>::quiet_NaN()};
  double vy{std::numeric_limits<double>::quiet_NaN()};
  double vz{std::numeric_limits<double>::quiet_NaN()};
  double t0{std::numeric_limits<double>::quiet_NaN()};

  //! id tag for this vertex
  int id{std::numeric_limits<int>::min()};

  ClassDefOverride(PHG4VtxPointv1, 2)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV1_H
#define G4DETECTORS_PHG4CYLINDERGEOMV1_H

#include "PHG4CylinderGeom.h"

#include <cmath>
#include <iostream>  // for cout, ostream

class PHParameters;

class PHG4CylinderGeomv1 : public PHG4CylinderGeom
{
 public:
  PHG4CylinderGeomv1() {}
  PHG4CylinderGeomv1(const double r, const double zmi, const double zma, const double thickn)
    : radius(r)
    , zmin(zmi)
    , zmax(zma)
    , thickness(thickn)
  {
  }

  ~PHG4CylinderGeomv1() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int get_layer() const override { return layer; }
  double get_radius() const override { return radius; }
  double get_thickness() const override { return thickness; }
  double get_zmin() const override { return zmin; }
  double get_zmax() const override { return zmax; }

  void set_layer(const int i) override { layer = i; }
  void set_radius(const double r) override { radius = r; }
  void set_thickness(const double t) override { thickness = t; }
  void set_zmin(const double z) override { zmin = z; }
  void set_zmax(const double z) override { zmax = z; }

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  void ImportParameters(const PHParameters& param) override;

 protected:
  int layer = -1;
  double radius = NAN;
  double zmin = NAN;
  double zmax = NAN;
  double thickness = NAN;

  ClassDefOverride(PHG4CylinderGeomv1, 1)
};

#endif

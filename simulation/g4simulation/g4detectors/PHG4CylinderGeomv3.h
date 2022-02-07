// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV3_H
#define G4DETECTORS_PHG4CYLINDERGEOMV3_H

#include "PHG4CylinderGeomv2.h"

#include <iostream>  // for cout, ostream

class PHG4CylinderGeomv3 : public PHG4CylinderGeomv2
{
 public:
  PHG4CylinderGeomv3();
  PHG4CylinderGeomv3(const double r, const double zmi, const double zma, const double thickn, const int n_scint,
                     const double tangl, const double phi_slat_null)
    : PHG4CylinderGeomv2(r, zmi, zma, thickn, n_scint)
    , tiltangle(tangl)
    , phi_slat_zero(phi_slat_null)
  {
  }

  ~PHG4CylinderGeomv3() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_tiltangle(const double phi) override { tiltangle = phi; }
  void set_phi_slat_zero(const double phi) override { phi_slat_zero = phi; }

  double get_phi_slat_zero() const override { return phi_slat_zero; }
  double get_tiltangle() const override { return tiltangle; }

 protected:
  double tiltangle;
  double phi_slat_zero;

  ClassDefOverride(PHG4CylinderGeomv3, 1)
};

#endif

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV2_H
#define G4DETECTORS_PHG4CYLINDERGEOMV2_H

#include "PHG4CylinderGeomv1.h"

#include <iostream>  // for cout, ostream

class PHParameters;

class PHG4CylinderGeomv2 : public PHG4CylinderGeomv1
{
 public:
  PHG4CylinderGeomv2() {}
  PHG4CylinderGeomv2(const double r, const double zmi, const double zma, const double thickn, const int n_scint)
    : PHG4CylinderGeomv1(r, zmi, zma, thickn)
    , nscint(n_scint)
  {
  }

  ~PHG4CylinderGeomv2() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_nscint(const int i) override { nscint = i; }
  int get_nscint() const override { return nscint; }

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  void ImportParameters(const PHParameters& param) override;

 protected:
  int nscint = -9999;

  ClassDefOverride(PHG4CylinderGeomv2, 1)
};

#endif

#ifndef PHFIELD_PHFIELD3DCARTESIAN_H
#define PHFIELD_PHFIELD3DCARTESIAN_H

#include "PHField.h"

#include <string>

//! untested code - I don't know if this is being used, drop me a line (with the field) and I test this - Chris P.
class PHField3DCartesian : public PHField
{
 public:
  PHField3DCartesian(const std::string &fname, const float magfield_rescale = 1.0);
  virtual ~PHField3DCartesian();

  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const;

 protected:
  std::string filename;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  double xstepsize;
  double ystepsize;
  double zstepsize;
  // these are updated in a const method
  // to cache previous values
  mutable double xyz[2][2][2][3];
  mutable double bf[2][2][2][3];
  mutable double xkey_save;
  mutable double ykey_save;
  mutable double zkey_save;
  mutable int cache_hits;
  mutable int cache_misses;
};

#endif

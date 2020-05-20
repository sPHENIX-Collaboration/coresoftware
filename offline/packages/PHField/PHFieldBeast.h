#ifndef PHFIELD_PHFIELDBEAST_H
#define PHFIELD_PHFIELDBEAST_H

#include "PHField.h"

#include <map>
#include <string>
#include <vector>

class BeastMagneticField;

class PHFieldBeast : public PHField
{
 public:
  PHFieldBeast(const std::string &filename, const int verb = 0, const float magfield_rescale = 1.0);
  virtual ~PHFieldBeast() {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const;

 private:
  BeastMagneticField *m_BeastMagneticField;
  float m_MagFieldScale;
};

#endif  // PHFIELD_PHFIELDBEAST_H

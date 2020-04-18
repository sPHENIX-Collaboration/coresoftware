#ifndef PHFIELD_PHFIELDCLEO_H
#define PHFIELD_PHFIELDCLEO_H

#include "PHField.h"

#include <map>
#include <string>
#include <vector>

class PHFieldCleo : public PHField
{
 public:
  PHFieldCleo(const std::string &filename, const int verb = 0, const float magfield_rescale = 1.0);
  virtual ~PHFieldCleo() {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const;

 private:
  std::vector<std::vector<std::vector<double> > > xField;
  std::vector<std::vector<std::vector<double> > > yField;
  std::vector<std::vector<std::vector<double> > > zField;

  // The dimensions of the table
  int nx, ny, nz;
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;

  float m_MagFieldScale;
};

#endif  // PHFIELD_PHFIELDCLEO_H

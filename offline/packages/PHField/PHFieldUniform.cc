#include "PHFieldUniform.h"

#include <Geant4/G4SystemOfUnits.hh>

PHFieldUniform::PHFieldUniform(
    double field_mag_x,
    double field_mag_y,
    double field_mag_z)
  : field_mag_x_(field_mag_x * tesla)
  , field_mag_y_(field_mag_y * tesla)
  , field_mag_z_(field_mag_z * tesla)
{
}

void PHFieldUniform::GetFieldValue(const double /*point*/ [4], double *Bfield) const
{
  Bfield[0] = field_mag_x_;
  Bfield[1] = field_mag_y_;
  Bfield[2] = field_mag_z_;

  return;
}

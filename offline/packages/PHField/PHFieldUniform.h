
#ifndef __PHFieldUniform_H__
#define __PHFieldUniform_H__

#include "PHField.h"

#include <map>
#include <string>
#include <vector>

class PHFieldUniform : public PHField
{
 public:
  //! construct field map in constant in units of Tesla
  PHFieldUniform(
      double field_mag_x,
      double field_mag_y,
      double field_mag_z);
  virtual ~PHFieldUniform() {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const;

  double get_field_mag_x() const
  {
    return field_mag_x_;
  }

  void set_field_mag_x(double fieldMagX)
  {
    field_mag_x_ = fieldMagX;
  }

  double get_field_mag_y() const
  {
    return field_mag_y_;
  }

  void set_field_mag_y(double fieldMagY)
  {
    field_mag_y_ = fieldMagY;
  }

  double get_field_mag_z() const
  {
    return field_mag_z_;
  }

  void set_field_mag_z(double fieldMagZ)
  {
    field_mag_z_ = fieldMagZ;
  }

 protected:
  double field_mag_x_;
  double field_mag_y_;
  double field_mag_z_;

 private:
};

#endif  // __PHFieldUniform_H

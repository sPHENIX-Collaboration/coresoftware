

#ifndef __PHFIELD_H__
#define __PHFIELD_H__

#include <phool/PHObject.h>

#include <vector>

//! \brief transient DST object for field storage and access
class PHField : public PHObject
{
 public:
  //! constructor
  explicit PHField(const int verb = 0)
    : verb_(verb)
  {
  }
  virtual ~PHField() {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in units of cm, s
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in units of Tesla
  virtual void GetFieldValue(
      const double Point[4],
      double *Bfield) const = 0;

 protected:
  void Verbosity(const int i) { verb_ = i; }
};

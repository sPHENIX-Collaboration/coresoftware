#ifndef PHFIELD_PHFIELD3DCARTESIAN_H
#define PHFIELD_PHFIELD3DCARTESIAN_H

#include "PHField.h"

#include <string>

//! untested code - I don't know if this is being used, drop me a line (with the field) and I test this - Chris P.
class PHField3DCartesian : public PHField
{
 public:
  PHField3DCartesian(const std::string &fname, const float magfield_rescale = 1.0);
  ~PHField3DCartesian() override;

  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const override;

template<typename T>
struct Interpolator
{
  /**
   * Linear interpolation at point p
   *
   *        (vx0)            (vx1)
   *     --o---------o------o------
   *       |         |      |      x
   *       0         p      1
   *
   */
  static T Linear(const T px, const T vx[2])
  {
    return vx[0] + (vx[1] - vx[0]) * px;
  }

  /**
   * Bilinear interpolation at point p
   *
   *     y
   *       |
   *       |
   *    1--o----------------o
   *       |(vy0)    |      |(vy1)
   *       |         |      |
   *   py--|---------p------|
   *       |         |      |
   *       |         |      |
   *       |         |      |
   *    0--o----------------o------
   *       |(vx0)    |      |(vx1) x
   *       0         px     1
   *
   */
  static T Bilinear(const T px, const T py, const T vx[2], const T vy[2])
  {
    T v_tmp[2] = { Linear(px, vx), Linear(px, vy), };
    return Linear(py, v_tmp);
  }

  /**
   * Trilinear interpolation at point p
   *
   *            o----------------o
   *           /|(vz2)          /|(vz3)
   *     y    / |              / |
   *       | /  |             /  |
   *       |/   |            /   |
   *    1--o----------------o    |
   *  (vy0)|    |/     (vy1)|    |
   *       | 1--o-----------|----o
   *       |   / (vz0)      |   / (vz1)
   *       |  /             |  /
   *       | /              | /
   *       |/               |/
   *    0--o----------------o------
   *      /|(vx0)           |(vx1) x
   *   z / 0                1
   *
   */
  static T Trilinear(const T px, const T py, const T pz, const T vx[2], const T vy[2], const T vz[4])
  {
    T v_tmp[2] = { Bilinear(px, py, vx, vy), Bilinear(px, py, vz, vz+2) };
    return Linear(pz, v_tmp);
  }
};


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

#ifndef PHFIELD_PHFIELD3DCARTESIAN_H
#define PHFIELD_PHFIELD3DCARTESIAN_H

#include "PHField.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <cmath>
#include <map>
#include <set>
#include <string>

class PHField3DCartesian : public PHField
{
 public:
explicit PHField3DCartesian(const std::string &fname, const float magfield_rescale = 1.0, const float innerradius = 0, const float outerradius = 1.e10, const float size_z = 1.e10);
  ~PHField3DCartesian() override;

  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const override;

 private:
  std::string filename;
  double xmin = 1000000;
  double xmax = -1000000;
  double ymin = 1000000;
  double ymax = -1000000;
  double zmin = 1000000;
  double zmax = -1000000;
  double xstepsize = NAN;
  double ystepsize = NAN;
  double zstepsize = NAN;
  // these are updated in a const method
  // to cache previous values
  mutable double xyz[2][2][2][3]{};
  mutable double bf[2][2][2][3]{};
  mutable double xkey_save = NAN;
  mutable double ykey_save = NAN;
  mutable double zkey_save = NAN;
  mutable int cache_hits = 0;
  mutable int cache_misses = 0;

  typedef boost::tuple<float, float, float> trio;
  std::map<boost::tuple<float, float, float>, boost::tuple<float, float, float> > fieldmap;
  std::set<float> xvals;
  std::set<float> yvals;
  std::set<float> zvals;
};

#endif

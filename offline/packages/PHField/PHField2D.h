
#ifndef PHFIELD_PHFIELD2D_H
#define PHFIELD_PHFIELD2D_H

#include "PHField.h"

#include <boost/tuple/tuple.hpp>

#include <map>
#include <string>
#include <vector>

class PHField2D : public PHField
{
  typedef boost::tuple<float, float> trio;

 public:
  PHField2D(const std::string &filename, const int verb = 0, const float magfield_rescale = 1.0);
  ~PHField2D() override {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  void GetFieldValue(const double Point[4], double *Bfield) const override;

  void GetFieldCyl(const double CylPoint[4], double *Bfield) const;

 protected:
  // < i, j, k > , this allows i and i+1 to be neighbors ( <i,j,k>=<z,r,phi> )
  std::vector<std::vector<float> > BFieldZ_;
  std::vector<std::vector<float> > BFieldR_;
  std::vector<std::vector<float> > BFieldPHI_;

  // maps indices to values z_map[i] = z_value that corresponds to ith index
  std::vector<float> z_map_;    // < i >
  std::vector<float> r_map_;    // < j >
  std::vector<float> phi_map_;  // < k >

  float maxz_, minz_;  // boundaries of magnetic field map cyl
  double magfield_unit;

 private:
  void print_map(std::map<trio, trio>::iterator &it) const;
  // mutable allows to change internal data even in const methods
  // I don't like this too much but these are cached values to speed up
  // the field lookup by a lot
  // I want them to be data members so we can run 2 fieldmaps in parallel
  // and still have caching. Putting those as static variables into
  // the implementation will prevent this
  mutable unsigned int r_index0_cache;
  mutable unsigned int r_index1_cache;
  mutable unsigned int z_index0_cache;
  mutable unsigned int z_index1_cache;
};

#endif

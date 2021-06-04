//
//    *************************************
//    *                                   *
//    *          PHField3DCylindrical.h            *
//    *                                   *
//    *************************************
//
// **************************************************************
//
// GEANT Field Map, for use in converting ROOT maps of
// PHENIX magnetic field.  Written by Michael Stone, July 2011
//
// **************************************************************
//
// The structure of the indices has little to do with physics and
// has much more to do with the way the PHENIX field map is formatted
// in SimMap3D++.root i.e.  The z value is incremented only after
// every phi and r point has been accounted for in that plane.

#ifndef PHFIELD_PHFIELD3DCYLINDRICAL_H
#define PHFIELD_PHFIELD3DCYLINDRICAL_H

#include "PHField.h"

#include <boost/tuple/tuple.hpp>

#include <map>
#include <string>
#include <vector>

class PHField3DCylindrical : public PHField
{
  typedef boost::tuple<float, float, float> trio;

 public:
  PHField3DCylindrical(const std::string& filename, int verb = 0, const float magfield_rescale = 1.0);
  ~PHField3DCylindrical() override {}
  void GetFieldValue(const double Point[4], double* Bfield) const override;
  void GetFieldCyl(const double CylPoint[4], double* Bfield) const;

 protected:
  // < i, j, k > , this allows i and i+1 to be neighbors ( <i,j,k>=<z,r,phi> )
  std::vector<std::vector<std::vector<float> > > BFieldZ_;
  std::vector<std::vector<std::vector<float> > > BFieldR_;
  std::vector<std::vector<std::vector<float> > > BFieldPHI_;

  // maps indices to values z_map[i] = z_value that corresponds to ith index
  std::vector<float> z_map_;    // < i >
  std::vector<float> r_map_;    // < j >
  std::vector<float> phi_map_;  // < k >

  float maxz_, minz_;  // boundaries of magnetic field map cyl

 private:
  bool bin_search(const std::vector<float>& vec, unsigned start, unsigned end, const float& key, unsigned& index) const;
  void print_map(std::map<trio, trio>::iterator& it) const;
};

#endif

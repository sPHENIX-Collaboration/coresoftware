/*
   Authors: Haiwang Yu
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_Field2D_h
#define genfit_Field2D_h

#include "GenFit/AbsBField.h"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include <map>

#include <TH2D.h>


namespace genfit {

/** @brief 2D Field
 *
 *  @author Haiwang Yu (New Mexico State University)
 * 
 */
class Field2D : public AbsBField {
 public:
  //! define the constant field in this ctor
  Field2D() : field_map_r_(NULL), field_map_z_(NULL)
  { ; }

  Field2D(std::string inname)
  { initialize(inname); }

  bool initialize(std::string inname = "");

  bool re_scale(double r);

  void plot(std::string option = "");

  //! return value at position
  TVector3 get(const TVector3& pos) const;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const;

 private:
  TH2D *field_map_r_;
  TH2D *field_map_z_;
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_Field2D_h

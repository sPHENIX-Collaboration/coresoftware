/*
   Authors: Haiwang Yu
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_Field_h
#define genfit_Field_h

#include "GenFit/AbsBField.h"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include <map>

#include <TH2D.h>

class PHField;

namespace genfit {

/** @brief Field Wrapper
 *
 *  @author Haiwang Yu (New Mexico State University)
 * 
 */
class Field : public AbsBField {
 public:
  Field(const PHField* field);

  virtual ~Field() {}


//  void plot(std::string option = "");

  //! return value at position
  TVector3 get(const TVector3& pos) const;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const;

  const PHField* get_field() const
  {
    return field_;
  }

  void set_field(const PHField* field)
  {
    field_ = field;
  }

 private:
  const PHField* field_;
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_Field_h

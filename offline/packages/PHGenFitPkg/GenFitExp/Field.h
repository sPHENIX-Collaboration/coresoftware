/*
   Authors: Haiwang Yu
*/
/** @addtogroup genfit
 * @{
 */

#ifndef GENFITEXP_FIELD_H
#define GENFITEXP_FIELD_H

#include <GenFit/AbsBField.h>

#include <TVector3.h>

class PHField;

namespace genfit
{
/** @brief Field Wrapper
 *
 *  @author Haiwang Yu (New Mexico State University)
 * 
 */
class Field : public AbsBField
{
 public:
  Field(const PHField* field);

  ~Field() override {}

  //  void plot(std::string option = "");

  //! return value at position
  TVector3 get(const TVector3& pos) const override;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const override;

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

#endif  // genfit_Field_h

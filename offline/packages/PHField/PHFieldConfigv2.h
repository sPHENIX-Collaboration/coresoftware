// $Id: $

/*!
 * \file PHFieldConfigv2.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFIELD_PHFIELDCONFIGV2_H
#define PHFIELD_PHFIELDCONFIGV2_H

#include "PHFieldConfig.h"

#include <iostream>

class PHObject;

/*!
 * \brief PHFieldConfigv2 implements field configuration information for uniform field model*/
class PHFieldConfigv2 : public PHFieldConfig
{
 public:
  //! construct field configuration in units of Tesla
  PHFieldConfigv2(
      double field_mag_x,
      double field_mag_y,
      double field_mag_z);

  //! default constructor for ROOT file IO
  PHFieldConfigv2()
    : PHFieldConfigv2(0, 0, 0)
  {
  }

  ~PHFieldConfigv2() override {}

  /// Virtual copy constructor.
  PHObject* CloneMe() const override { return new PHFieldConfigv2(*this); }

  /** identify Function from PHObject
   @param os Output Stream
   */
  void
  identify(std::ostream& os = std::cout) const override;

  /// Clear Content
  void Reset() override {}

  /// isValid returns non zero if object contains vailid data
  int
  isValid() const override { return 3; }

  FieldConfigTypes get_field_config() const override
  {
    return kFieldUniform;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_x() const override
  {
    return field_mag_x_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_x(double fieldMagX) override
  {
    field_mag_x_ = fieldMagX;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_y() const override
  {
    return field_mag_y_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_y(double fieldMagY) override
  {
    field_mag_y_ = fieldMagY;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_z() const override
  {
    return field_mag_z_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_z(double fieldMagZ) override
  {
    field_mag_z_ = fieldMagZ;
  }

 protected:
  double field_mag_x_;
  double field_mag_y_;
  double field_mag_z_;

  ClassDefOverride(PHFieldConfigv2, 1)
};

#endif

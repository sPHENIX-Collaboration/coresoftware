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

  virtual ~PHFieldConfigv2() {}

  /// Virtual copy constructor.
  virtual PHObject* CloneMe() const { return new PHFieldConfigv2(*this); }

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const;

  /// Clear Content
  virtual void
  Reset() {}

  /// isValid returns non zero if object contains vailid data
  virtual int
  isValid() const { return 3; }

  FieldConfigTypes get_field_config() const
  {
    return kFieldUniform;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_x() const
  {
    return field_mag_x_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_x(double fieldMagX)
  {
    field_mag_x_ = fieldMagX;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_y() const
  {
    return field_mag_y_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_y(double fieldMagY)
  {
    field_mag_y_ = fieldMagY;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  double get_field_mag_z() const
  {
    return field_mag_z_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfigv2
  void set_field_mag_z(double fieldMagZ)
  {
    field_mag_z_ = fieldMagZ;
  }

 protected:
  double field_mag_x_;
  double field_mag_y_;
  double field_mag_z_;

  ClassDef(PHFieldConfigv2, 1)
};

#endif

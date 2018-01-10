// $Id: $

/*!
 * \file PHFieldConfig_v2.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFieldConfig_v2_H_
#define PHFieldConfig_v2_H_

#include "PHFieldConfig.h"


/*!
 * \brief PHFieldConfig_v2 implements field configuration information for uniform field model*/
class PHFieldConfig_v2 : public PHFieldConfig
{
 public:
  //! construct field configuration in units of Tesla
  PHFieldConfig_v2(
      double field_mag_x,
      double field_mag_y,
      double field_mag_z
      );

  //! default constructor for ROOT file IO
  PHFieldConfig_v2(): PHFieldConfig_v2(0,0,0) {}

  virtual ~PHFieldConfig_v2();

  /// Virtual copy constructor.
  virtual PHObject*
  clone() const;

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void
  Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int
  isValid() const;

  FieldConfigTypes get_field_config() const
  {
    return kFieldUniform;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  double get_field_mag_x() const
  {
    return field_mag_x_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  void set_field_mag_x(double fieldMagX)
  {
    field_mag_x_ = fieldMagX;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  double get_field_mag_y() const
  {
    return field_mag_y_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  void set_field_mag_y(double fieldMagY)
  {
    field_mag_y_ = fieldMagY;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  double get_field_mag_z() const
  {
    return field_mag_z_;
  }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  void set_field_mag_z(double fieldMagZ)
  {
    field_mag_z_ = fieldMagZ;
  }

 protected:

  double field_mag_x_;
  double field_mag_y_;
  double field_mag_z_;

  ClassDef(PHFieldConfig_v2, 1)
};

#endif /* PHFieldConfig_v2_H_ */

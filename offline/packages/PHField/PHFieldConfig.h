// $Id: $

/*!
 * \file PHFieldConfig.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFIELD_PHFIELDCONFIG_H
#define PHFIELD_PHFIELDCONFIG_H

#include <phool/PHObject.h>

#include <iostream>
#include <string>

/*!
 * \brief PHFieldConfig store field configuration information */
class PHFieldConfig : public PHObject
{
 public:
  virtual ~PHFieldConfig(){}

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

  virtual PHObject* CloneMe() const {return nullptr;}

  enum FieldConfigTypes
  {
    //! Constant field
    kFieldUniform = 0,
    //! 2D field map expressed in cylindrical coordinates
    kField2D = 2,
    //! 3D field map expressed in cylindrical coordinates
    kField3DCylindrical = 3,
    //! 3D field map expressed in Cartesian coordinates
    Field3DCartesian = 1,

    //! invalid value
    kFieldInvalid = 9999
  };

  virtual FieldConfigTypes get_field_config() const;

  std::string get_field_config_description() const;

  virtual void set_field_config(FieldConfigTypes fieldConfig);

  virtual const std::string& get_filename() const;

  virtual void set_filename(const std::string& filename);

  virtual double get_magfield_rescale() const;

  virtual void set_magfield_rescale(double magfieldRescale);

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_x() const;

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_x(double fieldMagX);

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_y() const;

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_y(double fieldMagY);

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_z() const;
  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_z(double fieldMagZ);

 protected:
  //! pure virtual interface class. not for direct use
  PHFieldConfig(){}

  static const std::string kInvalid_FileName;

  ClassDef(PHFieldConfig, 1)
};

#endif

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
#include <limits>
#include <string>

/*!
 * \brief PHFieldConfig store field configuration information */
class PHFieldConfig : public PHObject
{
 public:
  ~PHFieldConfig() override {}

  /** identify Function from PHObject
   @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  enum FieldConfigTypes
  {
    //! Constant field
    kFieldUniform = 0,
    //! 2D field map expressed in cylindrical coordinates
    kField2D = 2,
    //! 3D field map expressed in cylindrical coordinates
    kField3DCylindrical = 3,
    //! Beast field map from https://github.com/eic/BeastMagneticField
    kFieldBeast = 4,
    //! Cleo field map from https://gitlab.com/eic/escalate/g4e/-/blob/master/SolenoidMag3D.TABLE
    kFieldCleo = 5,
    //! 3D field map expressed in Cartesian coordinates
    Field3DCartesian = 1,

    //! invalid value
    kFieldInvalid = 9999
  };

  virtual FieldConfigTypes get_field_config() const { return kFieldInvalid; }

  std::string get_field_config_description() const;

  virtual void set_field_config(FieldConfigTypes /*fieldConfig*/) { return; }

  virtual const std::string& get_filename() const { return kInvalid_FileName; }

  virtual void set_filename(const std::string& /*filename*/) { return; }

  virtual double get_magfield_rescale() const { return std::numeric_limits<double>::signaling_NaN(); }

  virtual void set_magfield_rescale(double /*magfieldRescale*/) { return; }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_x() const { return std::numeric_limits<double>::signaling_NaN(); }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_x(double /*fieldMagX*/) { return; }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_y() const { return std::numeric_limits<double>::signaling_NaN(); }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_y(double /*fieldMagY*/) { return; }

  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual double get_field_mag_z() const { return std::numeric_limits<double>::signaling_NaN(); }
  //! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
  virtual void set_field_mag_z(double /*fieldMagZ*/) { return; }

 protected:
  //! pure virtual interface class. not for direct use
  PHFieldConfig() {}

  static const std::string kInvalid_FileName;

  ClassDefOverride(PHFieldConfig, 1)
};

#endif

// $Id: $

/*!
 * \file PHFieldConfig.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfig.h"

#include <iostream>
#include <limits>

using namespace std;

const std::string PHFieldConfig::kInvalid_FileName("INVALID FILE");

string PHFieldConfig::get_field_config_description() const
{
  switch (get_field_config())
  {
  case kFieldUniform:
    return "Uniform field";
    break;
  case kField2D:
    return "2D field map expressed in cylindrical coordinates";
    break;
  case kField3DCylindrical:
    return "3D field map expressed in cylindrical coordinates";
    break;
  case Field3DCartesian:
    return "3D field map expressed in Cartesian coordinates";
    break;
  default:
    return "Invalid Field";
  }
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfig::identify(std::ostream& os) const
{
  os << "PHFieldConfig::identify - isValid() = " << isValid() << endl;
}

/// Clear Event
void PHFieldConfig::Reset()
{
}

/// isValid returns non zero if object contains vailid data
int PHFieldConfig::isValid() const
{
  return 0;
}

PHFieldConfig::FieldConfigTypes PHFieldConfig::get_field_config() const
{
  return kFieldInvalid;
}

void PHFieldConfig::set_field_config(PHFieldConfig::FieldConfigTypes fieldConfig)
{
}

const std::string& PHFieldConfig::get_filename() const
{
  return kInvalid_FileName;
}

void PHFieldConfig::set_filename(const std::string& filename)
{
}

double PHFieldConfig::get_magfield_rescale() const
{
  return std::numeric_limits<double>::signaling_NaN();
}

void PHFieldConfig::set_magfield_rescale(double magfieldRescale)
{
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
double PHFieldConfig::get_field_mag_x() const
{
  return std::numeric_limits<double>::signaling_NaN();
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
void PHFieldConfig::set_field_mag_x(double fieldMagX)
{
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
double PHFieldConfig::get_field_mag_y() const
{
  return std::numeric_limits<double>::signaling_NaN();
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
void PHFieldConfig::set_field_mag_y(double fieldMagY)
{
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
double PHFieldConfig::get_field_mag_z() const
{
  return std::numeric_limits<double>::signaling_NaN();
}

//! field value in Tesla for uniform field model ONLY for PHFieldConfig_v2
void PHFieldConfig::set_field_mag_z(double fieldMagZ)
{
}

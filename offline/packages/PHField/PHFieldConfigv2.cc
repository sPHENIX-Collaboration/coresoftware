// $Id: $

/*!
 * \file PHFieldConfigv2.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfigv2.h"

#include <iostream>
#include <string>

PHFieldConfigv2::PHFieldConfigv2(
    double field_mag_x,
    double field_mag_y,
    double field_mag_z)
  : field_mag_x_(field_mag_x)
  , field_mag_y_(field_mag_y)
  , field_mag_z_(field_mag_z)
{
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfigv2::identify(std::ostream& os) const
{
  os << "PHFieldConfigv2::identify -";
  if (isValid())
  {
    os << " Field type of [" << get_field_config_description();
    os << "] with field vector of ";
    os << "[" << get_field_mag_x();
    os << ", " << get_field_mag_y();
    os << ", " << get_field_mag_z();
    os << "] tesla";
  }
  else
  {
    os << "Empty";
  }
  os << std::endl;
}

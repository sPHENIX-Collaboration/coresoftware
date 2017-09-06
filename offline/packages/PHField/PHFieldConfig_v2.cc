// $Id: $

/*!
 * \file PHFieldConfig_v2.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfig_v2.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TMemFile.h>

#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

PHFieldConfig_v2::PHFieldConfig_v2(
    double field_mag_x,
    double field_mag_y,
    double field_mag_z)
  : field_mag_x_(field_mag_x )
  , field_mag_y_(field_mag_y )
  , field_mag_z_(field_mag_z)
{
}

PHFieldConfig_v2::~PHFieldConfig_v2()
{
}

/// Virtual copy constructor.
PHObject*
PHFieldConfig_v2::clone() const
{
  return new PHFieldConfig_v2(*this);
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfig_v2::identify(std::ostream& os) const
{
  os << "PHFieldConfig_v2::identify -";
  if (isValid())
  {
    os << " Field type of [" << get_field_config_description();
    os << "] with field vector of ";
    os <<"["<<get_field_mag_x();
    os <<", "<<get_field_mag_y();
    os <<", "<<get_field_mag_z();
    os <<"] tesla";
  }
  else
    os << "Empty";
  os << endl;
}
/// Clear Event
void PHFieldConfig_v2::Reset()
{
}

/// isValid returns non zero if object contains vailid data
int PHFieldConfig_v2::isValid() const
{
  return 3;
}

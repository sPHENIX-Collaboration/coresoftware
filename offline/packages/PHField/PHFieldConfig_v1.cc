// $Id: $

/*!
 * \file PHFieldConfig_v1.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfig_v1.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TMemFile.h>

#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

PHFieldConfig_v1::PHFieldConfig_v1(FieldConfigTypes field_config,
                                   const std::string& filename,
                                   double magfield_rescale)
  : field_config_(field_config)
  , filename_(filename)
  , magfield_rescale_(magfield_rescale)
{
}

PHFieldConfig_v1::~PHFieldConfig_v1()
{
}

/// Virtual copy constructor.
PHObject*
PHFieldConfig_v1::clone() const
{
  return new PHFieldConfig_v1(*this);
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfig_v1::identify(std::ostream& os) const
{
  os << "PHFieldConfig_v1::identify -";
  if (isValid())
  {
    os << " Field type of [" << get_field_config_description();
    os << "] from file [" << get_filename();
    os << "] with a scale factor of " << get_magfield_rescale();
  }
  else
    os << "Empty";
  os << endl;
}
/// Clear Event
void PHFieldConfig_v1::Reset()
{
}

/// isValid returns non zero if object contains vailid data
int PHFieldConfig_v1::isValid() const
{
  return filename_.length();
}

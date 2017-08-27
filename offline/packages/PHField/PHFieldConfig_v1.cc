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
                                   std::string& filename,
                                   double magfield_rescale)
  : field_config_(field_config)
  , filename_(filename)
  , magfield_rescale_(magfield_rescale)
{
}

PHFieldConfig_v1::~PHFieldConfig_v1()
{
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfig_v1::identify(std::ostream& os) const
{
  os << "PHFieldConfig_v1::identify - ";
  if (isValid())
  {
    os << "\tget_field_config() \t= " << get_field_config() << enld;
    os << "\tget_filename() \t= " << get_filename() << enld;
    os << "\tget_magfield_rescale() \t= " << get_magfield_rescale() << enld;
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

// $Id: $

/*!
 * \file PHFieldConfig.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfig.h"

#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

const std::string PHFieldConfig::kInvalid_FileName("INVALID FILE");

PHFieldConfig::PHFieldConfig()
  : Data(0)
{
}

PHFieldConfig::~PHFieldConfig()
{
}

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

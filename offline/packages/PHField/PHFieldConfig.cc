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

PHFieldConfig::PHFieldConfig()
  : Data(0)
{
}

PHFieldConfig::~PHFieldConfig()
{
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfig::identify(std::ostream& os) const
{
  os << "PHFieldConfig - ";
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
void PHFieldConfig::Reset()
{
}

/// isValid returns non zero if object contains vailid data
int PHFieldConfig::isValid() const
{
  return 0;
}

// $Id: $

/*!
 * \file PHFieldConfigv1.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHFieldConfigv1.h"

#include <iostream>

PHFieldConfigv1::PHFieldConfigv1(FieldConfigTypes field_config,
                                 const std::string& filename,
                                 double magfield_rescale)
  : field_config_(field_config)
  , filename_(filename)
  , magfield_rescale_(magfield_rescale)
{
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHFieldConfigv1::identify(std::ostream& os) const
{
  os << "PHFieldConfigv1::identify -";
  if (isValid())
  {
    os << " Field type of [" << get_field_config_description();
    os << "] from file [" << get_filename();
    os << "] with a scale factor of " << get_magfield_rescale();
  }
  else
  {
    os << "Empty";
  }
  os << std::endl;
}

/// isValid returns non zero if object contains valid data
int PHFieldConfigv1::isValid() const
{
  return filename_.length();
}

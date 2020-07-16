#include "EicEventHeaderv1.h"

#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, basic_ostream::o...
#include <string>    // for operator<<, string, char_traits
#include <utility>   // for pair

using namespace std;

EicEventHeaderv1::EicEventHeaderv1(const EicEventHeader *eicevt)
{
  CopyFrom(eicevt);
}

void EicEventHeaderv1::Reset()
{
  prop_map.clear();
}

bool EicEventHeaderv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

float EicEventHeaderv1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_float))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_float) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).fdata;

  return NAN;
}

int EicEventHeaderv1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_int))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_int) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).idata;

  return INT_MIN;
}

unsigned int
EicEventHeaderv1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_uint))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_uint) << endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).uidata;

  return UINT_MAX;
}

void EicEventHeaderv1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id, type_float))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_float) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void EicEventHeaderv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_int) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void EicEventHeaderv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    cout << PHWHERE << " Property " << property_info.first << " with id "
         << prop_id << " is of type " << get_property_type(property_info.second)
         << " not " << get_property_type(type_uint) << endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
EicEventHeaderv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
  {
    return iter->second;
  }
  return UINT_MAX;
}

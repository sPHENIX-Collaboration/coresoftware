#include "PHG4Cellv1.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <phool/phool.h>

#include <climits>  // for UINT_MAX, INT_MIN
#include <cmath>    // for NAN
#include <cstdlib>  // for exit
#include <iostream>
#include <string>  // for operator<<, string

PHG4Cellv1::PHG4Cellv1(const PHG4CellDefs::keytype g4cellid)
  : cellid(g4cellid)
{
}

PHG4Cellv1::~PHG4Cellv1()
{
  hitedeps.clear();
  showeredeps.clear();
  prop_map.clear();
  return;
}

void PHG4Cellv1::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep)
{
  hitedeps[g4hitid] = edep;
  return;
}

void PHG4Cellv1::add_shower_edep(const int g4showerid, const float edep)
{
  showeredeps[g4showerid] += edep;
  return;
}

bool PHG4Cellv1::has_binning(const PHG4CellDefs::CellBinning binning) const
{
  return PHG4CellDefs::has_binning(cellid, binning);
}

short int
PHG4Cellv1::get_detid() const
{
  return PHG4CellDefs::get_detid(cellid);
}

bool PHG4Cellv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

float PHG4Cellv1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_float))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_float) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).fdata;

  return NAN;
}

int PHG4Cellv1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).idata;

  return INT_MIN;
}

unsigned int
PHG4Cellv1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end()) return u_property(i->second).uidata;

  return UINT_MAX;
}

void PHG4Cellv1::add_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id, type_float))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_float) << std::endl;
    exit(1);
  }
  float val = value;
  if (prop_map.find(prop_id) != prop_map.end())
  {
    val += get_property_float(prop_id);
  }
  prop_map[prop_id] = u_property(val).get_storage();
}

void PHG4Cellv1::add_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  int val = value;
  if (prop_map.find(prop_id) != prop_map.end())
  {
    val += get_property_int(prop_id);
  }
  prop_map[prop_id] += u_property(val).get_storage();
}

void PHG4Cellv1::add_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  unsigned int val = value;
  if (prop_map.find(prop_id) != prop_map.end())
  {
    val += get_property_uint(prop_id);
  }
  prop_map[prop_id] += u_property(val).get_storage();
}

void PHG4Cellv1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id, type_float))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_float) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void PHG4Cellv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  prop_map[prop_id] += u_property(value).get_storage();
}

void PHG4Cellv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  prop_map[prop_id] += u_property(value).get_storage();
}

unsigned int
PHG4Cellv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
  {
    return iter->second;
  }
  return UINT_MAX;
}

void PHG4Cellv1::print() const
{
  identify(std::cout);
}

void PHG4Cellv1::Reset()
{
  hitedeps.clear();
  showeredeps.clear();
  prop_map.clear();
  return;
}

void PHG4Cellv1::identify(std::ostream& os) const
{
  os << "PHG4Cellv1  0x" << std::hex << cellid << std::dec << std::endl;

  os << "Associated to " << hitedeps.size() << " hits" << std::endl;
  for (const auto pair : hitedeps)
  {
    os << "\t PHG4hit " << pair.first << " -> " << pair.second << " GeV" << std::endl;
  }

  os << "Associated to " << showeredeps.size() << " showers" << std::endl;
  for (const auto pair : showeredeps)
  {
    os << "\t Shower " << pair.first << " -> " << pair.second << " GeV" << std::endl;
  }
  os << "Properties:" << std::endl;
  for (auto i : prop_map)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i.first);
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    os << "\t" << prop_id << ":\t" << property_info.first << " = \t";
    switch (property_info.second)
    {
    case type_int:
      os << get_property_int(prop_id);
      break;
    case type_uint:
      os << get_property_uint(prop_id);
      break;
    case type_float:
      os << get_property_float(prop_id);
      break;
    default:
      os << " unknown type ";
    }
    os << std::endl;
  }
}

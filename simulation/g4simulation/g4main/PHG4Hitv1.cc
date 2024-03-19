#include "PHG4Hitv1.h"
#include "PHG4HitDefs.h"

#include <phool/phool.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <utility>

PHG4Hitv1::PHG4Hitv1(const PHG4Hit* g4hit)
{
  CopyFrom(g4hit);
}

void PHG4Hitv1::Reset()
{
  hitid = std::numeric_limits<PHG4HitDefs::keytype>::max();
  trackid = std::numeric_limits<int>::min();
  showerid = std::numeric_limits<int>::min();
  edep = std::numeric_limits<float>::quiet_NaN();
  for (int i = 0; i < 2; i++)
  {
    set_x(i, std::numeric_limits<float>::quiet_NaN());
    set_y(i, std::numeric_limits<float>::quiet_NaN());
    set_z(i, std::numeric_limits<float>::quiet_NaN());
    set_t(i, std::numeric_limits<float>::quiet_NaN());
  }
  prop_map.clear();
}

int PHG4Hitv1::get_detid() const
{
  // a compile time check if the hit_idbits are within range (1-32)
  static_assert(PHG4HitDefs::hit_idbits <= sizeof(unsigned int) * 8, "hit_idbits < 32, fix in PHG4HitDefs.h");
  int detid = (hitid >> PHG4HitDefs::hit_idbits);
  return detid;
}

void PHG4Hitv1::print() const
{
  std::cout << "New Hitv1  0x" << std::hex << hitid
            << std::dec << "  on track " << trackid << " EDep " << edep << std::endl;
  std::cout << "Location: X " << x[0] << "/" << x[1] << "  Y " << y[0] << "/" << y[1] << "  Z " << z[0] << "/" << z[1] << std::endl;
  std::cout << "Time        " << t[0] << "/" << t[1] << std::endl;

  for (auto i : prop_map)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i.first);
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << "\t" << prop_id << ":\t" << property_info.first << " = \t";
    switch (property_info.second)
    {
    case type_int:
      std::cout << get_property_int(prop_id);
      break;
    case type_uint:
      std::cout << get_property_uint(prop_id);
      break;
    case type_float:
      std::cout << get_property_float(prop_id);
      break;
    default:
      std::cout << " unknown type ";
    }
    std::cout << std::endl;
  }
}

bool PHG4Hitv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

float PHG4Hitv1::get_property_float(const PROPERTY prop_id) const
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

  if (i != prop_map.end())
  {
    return u_property(i->second).fdata;
  }

  return std::numeric_limits<float>::quiet_NaN();
}

int PHG4Hitv1::get_property_int(const PROPERTY prop_id) const
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

  if (i != prop_map.end())
  {
    return u_property(i->second).idata;
  }

  return std::numeric_limits<int>::min();
}

unsigned int
PHG4Hitv1::get_property_uint(const PROPERTY prop_id) const
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

  if (i != prop_map.end())
  {
    return u_property(i->second).uidata;
  }

  return std::numeric_limits<unsigned int>::max();
}

void PHG4Hitv1::set_property(const PROPERTY prop_id, const float value)
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

void PHG4Hitv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void PHG4Hitv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
PHG4Hitv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
  {
    return iter->second;
  }
  return std::numeric_limits<unsigned int>::max();
}

float PHG4Hitv1::get_px(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_px_0);
  case 1:
    return get_property_float(prop_px_1);
  default:
    std::cout << "Invalid index in get_px: " << i << std::endl;
    exit(1);
  }
}

float PHG4Hitv1::get_py(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_py_0);
  case 1:
    return get_property_float(prop_py_1);
  default:
    std::cout << "Invalid index in get_py: " << i << std::endl;
    exit(1);
  }
}

float PHG4Hitv1::get_pz(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_pz_0);
  case 1:
    return get_property_float(prop_pz_1);
  default:
    std::cout << "Invalid index in get_pz: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_px(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_px_0, f);
    return;
  case 1:
    set_property(prop_px_1, f);
    return;
  default:
    std::cout << "Invalid index in set_px: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_py(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_py_0, f);
    return;
  case 1:
    set_property(prop_py_1, f);
    return;
  default:
    std::cout << "Invalid index in set_py: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_pz(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_pz_0, f);
    return;
  case 1:
    set_property(prop_pz_1, f);
    return;
  default:
    std::cout << "Invalid index in set_pz: " << i << std::endl;
    exit(1);
  }
}

float PHG4Hitv1::get_local_x(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_local_x_0);
  case 1:
    return get_property_float(prop_local_x_1);
  default:
    std::cout << "Invalid index in get_local_x: " << i << std::endl;
    exit(1);
  }
}

float PHG4Hitv1::get_local_y(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_local_y_0);
  case 1:
    return get_property_float(prop_local_y_1);
  default:
    std::cout << "Invalid index in get_local_y: " << i << std::endl;
    exit(1);
  }
}

float PHG4Hitv1::get_local_z(const int i) const
{
  switch (i)
  {
  case 0:
    return get_property_float(prop_local_z_0);
  case 1:
    return get_property_float(prop_local_z_1);
  default:
    std::cout << "Invalid index in get_local_z: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_local_x(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_local_x_0, f);
    return;
  case 1:
    set_property(prop_local_x_1, f);
    return;
  default:
    std::cout << "Invalid index in set_local_x: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_local_y(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_local_y_0, f);
    return;
  case 1:
    set_property(prop_local_y_1, f);
    return;
  default:
    std::cout << "Invalid index in set_local_y: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::set_local_z(const int i, const float f)
{
  switch (i)
  {
  case 0:
    set_property(prop_local_z_0, f);
    return;
  case 1:
    set_property(prop_local_z_1, f);
    return;
  default:
    std::cout << "Invalid index in set_local_z: " << i << std::endl;
    exit(1);
  }
}

void PHG4Hitv1::identify(std::ostream& os) const
{
  os << "Class " << this->ClassName() << std::endl;
  os << "hitid: 0x" << std::hex << hitid << std::dec << std::endl;
  os << "x0: " << get_x(0)
     << ", y0: " << get_y(0)
     << ", z0: " << get_z(0)
     << ", t0: " << get_t(0) << std::endl;
  os << "x1: " << get_x(1)
     << ", y1: " << get_y(1)
     << ", z1: " << get_z(1)
     << ", t1: " << get_t(1) << std::endl;
  os << "trackid: " << trackid << ", showerid: " << showerid
     << ", edep: " << edep << std::endl;
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

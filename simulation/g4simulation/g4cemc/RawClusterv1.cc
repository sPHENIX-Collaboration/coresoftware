#include "RawClusterv1.h"
#include <limits>

using namespace std;

RawClusterv1::RawClusterv1()
  : RawCluster()
  , clusterid(0)
  , _z(numeric_limits<float>::signaling_NaN())
  , _r(numeric_limits<float>::signaling_NaN())
  , _phi(numeric_limits<float>::signaling_NaN())
  , _energy(numeric_limits<float>::signaling_NaN())
{
}

void RawClusterv1::Reset()
{
  clusterid = 0;
  _z = (numeric_limits<float>::signaling_NaN());
  _r = (numeric_limits<float>::signaling_NaN());
  _phi = (numeric_limits<float>::signaling_NaN());
  _energy = (numeric_limits<float>::signaling_NaN());
  towermap.clear();
}

void RawClusterv1::addTower(const RawClusterDefs::keytype twrid, const float etower)
{
  if (towermap.find(twrid) != towermap.end())
  {
    cout << "tower 0x" << hex << twrid << ", dec: " << dec
         << twrid << " already exists, that is bad" << endl;
    exit(1);
  }
  towermap[twrid] = etower;
}

void RawClusterv1::identify(std::ostream& os = std::cout) const
{
  os << "RawClusterv1" << std::endl;
}

//! convert cluster location to psuedo-rapidity given a user chosen z-location
float RawClusterv1::get_eta(const float z) const
{
  if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
  return asinh((get_z() - z) / get_r());
}


//! convert cluster E_T given a user chosen z-location
float RawClusterv1::get_et(const float z) const
{
  if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
  return get_energy() * sin;
}

bool
RawClusterv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i!=prop_map.end();
}

float
RawClusterv1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_float) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).fdata;

  return   NAN ;
}

int
RawClusterv1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_int) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).idata;

  return INT_MIN;
}

unsigned int
RawClusterv1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_uint) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).uidata;

  return UINT_MAX ;
}

void
RawClusterv1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_float) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
RawClusterv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_int) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
RawClusterv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_uint) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
RawClusterv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
    {
      return iter->second;
    }
  return UINT_MAX;
}

#include "CentralityInfov1.h"

#include <limits>

void CentralityInfov1::identify(std::ostream& os) const
{
  os << "CentralityInfo: " << std::endl;

  return;
}

void CentralityInfov1::Reset()
{
  _quantity_map.clear();
  _centile_map.clear();
}

bool CentralityInfov1::has_quantity(const PROP prop_id) const
{
  return _quantity_map.find(prop_id) != _quantity_map.end();
}

void CentralityInfov1::set_quantity(const PROP prop_id, float value)
{
  _quantity_map[prop_id] = value;
}

float CentralityInfov1::get_quantity(const PROP prop_id) const
{
  if (!has_quantity(prop_id))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    return _quantity_map.at(prop_id);
  }
}

bool CentralityInfov1::has_centile(const PROP prop_id) const
{
  return _centile_map.find(prop_id) != _centile_map.end();
}

void CentralityInfov1::set_centile(const PROP prop_id, float value)
{
  _centile_map[prop_id] = value;
}

float CentralityInfov1::get_centile(const PROP prop_id) const
{
  if (!has_centile(prop_id))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    return _centile_map.at(prop_id);
  }
}

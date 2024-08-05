#include "RawTowerv2.h"

#include <iostream>
#include <limits>
#include <string>   // for operator<<, string
#include <utility>  // for pair

RawTowerv2::RawTowerv2(const RawTower& tower)
  : RawTowerv1(tower)
{
  // This is a generic copy of ALL properties a hit has
  // do not add explicit copies, they will be added to
  // the new hits with their default value increasing memory use
  for (unsigned int ic = 0; ic < std::numeric_limits<unsigned char>::max(); ic++)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(ic);
    if (tower.has_property(prop_id))
    {
      set_property(prop_id, tower.get_property(prop_id));
    }
  }
}

RawTowerv2::RawTowerv2(RawTowerDefs::keytype id)
  : RawTowerv1(id)
{
}

RawTowerv2::RawTowerv2(const unsigned int ieta, const unsigned int iphi)
  : RawTowerv1(ieta, iphi)
{
}

RawTowerv2::RawTowerv2(const RawTowerDefs::CalorimeterId caloid,
                       const unsigned int ieta, const unsigned int iphi)
  : RawTowerv1(caloid, ieta, iphi)
{
}

void RawTowerv2::Reset()
{
  RawTowerv1::Reset();
  prop_map.clear();
}

int RawTowerv2::isValid() const
{
  return RawTowerv1::isValid();
}

void RawTowerv2::identify(std::ostream& os) const
{
  os << "RawTowerv2: etabin: " << get_bineta() << ", phibin: " << get_binphi()
     << " energy=" << get_energy() << std::endl;

  for (auto i : prop_map)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i.first);
    const std::string property_info = get_property_info(prop_id);
    std::cout << "\t" << prop_id << ":\t" << property_info << " = \t" << get_property(prop_id) << std::endl;
  }
}

bool RawTowerv2::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

double
RawTowerv2::get_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end())
  {
    return i->second;
  }

  return std::numeric_limits<double>::signaling_NaN();
}

void RawTowerv2::set_property(const PROPERTY prop_id, const double value)
{
  prop_map[prop_id] = value;
}

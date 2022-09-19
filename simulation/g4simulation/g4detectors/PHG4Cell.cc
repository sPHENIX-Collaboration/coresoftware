#include "PHG4Cell.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <phool/PHObject.h>  // for PHObject

#include <cassert>
#include <cstdlib>

void PHG4Cell::CopyFrom(const PHObject* phobj)
{
  const PHG4Cell* g4cell = dynamic_cast<const PHG4Cell*>(phobj);
  assert(g4cell);
  set_cellid(g4cell->get_cellid());
  for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
  {
    PROPERTY prop_id = static_cast<PHG4Cell::PROPERTY>(ic);
    if (g4cell->has_property(prop_id))
    {
      set_property_nocheck(prop_id, g4cell->get_property_nocheck(prop_id));
    }
  }
}

void PHG4Cell::identify(std::ostream& os) const
{
  os << "Class " << this->ClassName() << std::endl;
  return;
}

std::ostream& operator<<(std::ostream& stream, const PHG4Cell* /*cell*/)
{
  stream << "PHG4Cell" << std::endl;
  return stream;
}

PHG4Cell::EdepConstRange PHG4Cell::get_g4hits()
{
  static std::map<PHG4HitDefs::keytype, float> dummy;
  return std::make_pair(dummy.begin(), dummy.end());
}

PHG4Cell::ShowerEdepConstRange PHG4Cell::get_g4showers()
{
  static std::map<int, float> dummy;
  return std::make_pair(dummy.begin(), dummy.end());
}

void PHG4Cell::Reset()
{
  std::cout << "Reset not implemented by daughter class" << std::endl;
  return;
}

std::pair<const std::string, PHG4Cell::PROPERTY_TYPE>
PHG4Cell::get_property_info(const PROPERTY prop_id)
{
  switch (prop_id)
  {
  case prop_stave_index:
    return std::make_pair("stave index", PHG4Cell::type_int);
  case prop_half_stave_index:
    return std::make_pair("half stave index", PHG4Cell::type_int);
  case prop_module_index:
    return std::make_pair("module index", PHG4Cell::type_int);
  case prop_chip_index:
    return std::make_pair("chip index", PHG4Cell::type_int);
  case prop_pixel_index:
    return std::make_pair("pixel index", PHG4Cell::type_int);
  case prop_phibin:
    return std::make_pair("phibin", PHG4Cell::type_int);
  case prop_zbin:
    return std::make_pair("zbin", PHG4Cell::type_int);
  case prop_ladder_z_index:
    return std::make_pair("ladder z index", PHG4Cell::type_int);
  case prop_ladder_phi_index:
    return std::make_pair("ladder phi index", PHG4Cell::type_int);
  case prop_edep:
    return std::make_pair("energy deposition", PHG4Cell::type_float);
  case prop_eion:
    return std::make_pair("ionizing energy loss", PHG4Cell::type_float);
  case prop_light_yield:
    return std::make_pair("light yield", PHG4Cell::type_float);
  case prop_raw_light_yield:
    return std::make_pair("raw light yield", PHG4Cell::type_float);
  default:
    std::cout << "PHG4Cell::get_property_info - Fatal Error - unknown index " << prop_id << std::endl;
    exit(1);
  }
}

bool PHG4Cell::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
  {
    return false;
  }
  return true;
}

std::string PHG4Cell::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch (prop_type)
  {
  case type_int:
    return "int";
  case type_uint:
    return "unsigned int";
  case type_float:
    return "float";
  default:
    return "unkown";
  }
}

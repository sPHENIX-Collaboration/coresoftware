#include "TrackerHitv1.h"

#include <phool/phool.h>

#include <iostream>

TrackerHitv1::TrackerHitv1()
  : hitid(TrackerDefs::KEYMAX)
{
}

TrackerHitv1::TrackerHitv1(TrackerDefs::keytype id)
  : hitid(id)
{
}

TrackerHitv1::~TrackerHitv1()
{
  hitedeps.clear();
  // prop_map.clear();
  return;
}

void TrackerHitv1::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep)
{
  hitedeps[g4hitid] = edep;
  return;
}

void TrackerHitv1::print() const
{
  identify(std::cout);
}

void TrackerHitv1::Reset()
{
  hitedeps.clear();
  // prop_map.clear();
  return;
}

void TrackerHitv1::identify(std::ostream& os) const
{
  os << "New TrackerHitv1  0x" << std::hex << get_hitid() << std::dec << std::endl;

  os << "Associated to " << hitedeps.size() << " hits" << std::endl;
  for (const auto pair : hitedeps)
  {
    os << "\t PHG4hit " << pair.first << " -> " << pair.second << " GeV" << std::endl;
  }

  // for (prop_map_t::const_iterator i = prop_map.begin(); i != prop_map.end(); ++i)
  // {
  //   PROPERTY prop_id = static_cast<PROPERTY>(i->first);
  //   pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
  //   os << "\t" << prop_id << ":\t" << property_info.first << " = \t";
  //   switch (property_info.second)
  //   {
  //   case type_int:
  //     os << get_property_int(prop_id);
  //     break;
  //   case type_uint:
  //     os << get_property_uint(prop_id);
  //     break;
  //   case type_float:
  //     os << get_property_float(prop_id);
  //     break;
  //   default:
  //     os << " unknown type ";
  //   }
  //   os << endl;
  // }
}

#include "InttEventInfoContainerv1.h"

InttEventInfoContainerv1::InttEventInfoContainerv1()
{
  Reset();
}

InttEventInfoContainerv1::~InttEventInfoContainerv1()
{
  //Do nothing
}

void
InttEventInfoContainerv1::identify(std::ostream& os) const
{
  os << "InttEventInfoContainer::identify()" << std::endl;
  os << "\tVersion 1" << std::endl;
}

void
InttEventInfoContainerv1::Reset()
{
  info_map.clear();
}

void
InttEventInfoContainerv1::AddInfo(KEY_t const& key, VAL_t const& val)
{
  info_map[key] = val;
}

InttEventInfoContainerv1::VAL_t&
InttEventInfoContainerv1::GetInfo(KEY_t const& key)
{
  return info_map[key];
}

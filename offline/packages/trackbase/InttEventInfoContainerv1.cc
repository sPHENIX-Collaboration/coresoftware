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
  os << "InttEventInfoContainer::identify" << std::endl;
  os << "\tVersion 1" << std::endl;
}

void
InttEventInfoContainerv1::Reset()
{
  bco_full = 0;
}

uint64_t
InttEventInfoContainerv1::get_bco_full() const
{
  return bco_full;
}

void
InttEventInfoContainerv1::set_bco_full(uint64_t const& _bco_full)
{
  bco_full = _bco_full;
}

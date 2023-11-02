#include "InttEventInfov1.h"

InttEventInfov1::InttEventInfov1()
{
  bco_full = 0;
}

InttEventInfov1::~InttEventInfov1()
{
  //Do nothing
}

void
InttEventInfov1::identify(std::ostream& os) const
{
  os << "InttEventInfo::identify" << std::endl;
  os << "\tVersion 1" << std::endl;
}

void
InttEventInfov1::Reset()
{
  bco_full = 0;
}

uint64_t
InttEventInfov1::get_bco_full() const
{
  return bco_full;
}

void
InttEventInfov1::set_bco_full(uint64_t const& _bco_full)
{
  bco_full = _bco_full;
}

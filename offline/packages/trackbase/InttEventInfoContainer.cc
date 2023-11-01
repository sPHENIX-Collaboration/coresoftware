#include "InttEventInfoContainer.h"

InttEventInfoContainer::InttEventInfoContainer()
{
  //Do nothing
}

InttEventInfoContainer::~InttEventInfoContainer()
{
  //Do nothing
}

void
InttEventInfoContainer::identify(std::ostream& os) const
{
  os << "InttEventInfoContainer::identify" << std::endl;
  os << "\tBase instance" << std::endl;
}

void
InttEventInfoContainer::Reset()
{
  std::cout << "InttEventInfoContainer::Reset" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);
}

uint64_t
InttEventInfoContainer::get_bco_full() const
{
  std::cout << "InttEventInfoContainer::get_bco_full" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);

  return 0;
}

void
InttEventInfoContainer::set_bco_full(uint64_t const&)
{
  std::cout << "InttEventInfoContainer::set_bco_full" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);
}

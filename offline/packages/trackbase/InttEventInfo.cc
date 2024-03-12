#include "InttEventInfo.h"

void InttEventInfo::identify(std::ostream& os) const
{
  os << "InttEventInfo::identify" << std::endl;
  os << "\tBase instance" << std::endl;
}

void InttEventInfo::Reset()
{
  std::cout << "InttEventInfo::Reset" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);
}

uint64_t
InttEventInfo::get_bco_full() const
{
  std::cout << "InttEventInfo::get_bco_full" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);

  return 0;
}

void InttEventInfo::set_bco_full(uint64_t const& /*unused*/)
{
  std::cout << "InttEventInfo::set_bco_full" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);
}

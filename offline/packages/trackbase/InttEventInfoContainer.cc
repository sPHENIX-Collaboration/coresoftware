#include "InttEventInfoContainer.h"

namespace
{
  InttEventInfoContainer::Map dummy_map;
}

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
  os << "InttEventInfoContainer::identify()" << std::endl;
  os << "\tBase instance" << std::endl;
}

void
InttEventInfoContainer::Reset()
{
  std::cout << "InttEventInfoContainer::Reset()" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);
}

void
InttEventInfoContainer::AddInfo(KEY_t const& key, VAL_t const& val)
{
  std::cout << "InttEventInfoContainer::AddInfo(KEY_t const&, VAL_t const&)" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);

  dummy_map[key] = val;
}

InttEventInfoContainer::VAL_t&
InttEventInfoContainer::GetInfo(KEY_t const& key)
{
  std::cout << "InttEventInfoContainer::GetInfo(KEY_t const&, VAL_t const&)" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  exit(1);

  return dummy_map[key];
}

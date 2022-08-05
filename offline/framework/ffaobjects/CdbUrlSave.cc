#include "CdbUrlSave.h"

#include <phool/phool.h>

#include <cstdint>  // for uint64_t
#include <iostream>

class PHObject;

static std::vector<std::tuple<std::string, std::string, uint64_t>> dummy;

PHObject*
CdbUrlSave::CloneMe() const
{
  std::cout << "CdbUrlSave::CloneMe() is not implemented in daugther class" << std::endl;
  return nullptr;
}

void CdbUrlSave::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void CdbUrlSave::identify(std::ostream& os) const
{
  os << "identify yourself: virtual CdbUrlSave Object" << std::endl;
  return;
}

int CdbUrlSave::isValid() const
{
  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
  return 0;
}

std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator CdbUrlSave::begin() const
{
  return dummy.begin();
}

std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator CdbUrlSave::end() const
{
  return dummy.begin();
}

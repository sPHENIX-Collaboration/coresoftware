#include "BcoInfo.h"

#include <phool/phool.h>

#include <iostream>

void BcoInfo::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void BcoInfo::identify(std::ostream& os) const
{
  os << "identify yourself: virtual BcoInfo Object" << std::endl;
  return;
}

int BcoInfo::isValid() const
{
  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
  return 0;
}

#include "StreamingBcoInfo.h"

#include <phool/phool.h>

#include <iostream>

void StreamingBcoInfo::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void StreamingBcoInfo::identify(std::ostream& os) const
{
  os << "identify yourself: virtual StreamingBcoInfo Object" << std::endl;
  return;
}

//int BcoStreamingLumiInfo::isValid() const
//{
//  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
//  return 0;
//}

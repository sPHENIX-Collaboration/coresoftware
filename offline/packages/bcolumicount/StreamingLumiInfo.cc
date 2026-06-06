#include "StreamingLumiInfo.h"

#include <phool/phool.h>

#include <iostream>

void StreamingLumiInfo::identify(std::ostream& os) const
{
  os << "identify yourself: virtual StreamingLumiInfo Object" << std::endl;
  return;
}

//int StreamingLumiInfo::isValid() const
//{
//  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
//  return 0;
//}

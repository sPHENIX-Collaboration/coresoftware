#include "StreamingLumiInfov1.h"

#include <phool/phool.h>

#include <iostream>

void StreamingLumiInfov1::identify(std::ostream& os) const
{
  os << "identify yourself: I am a StreamingLumiInfov1 Object\n";  return;
}

//int StreamingLumiInfov1::isValid() const
//{
//  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
//  return 0;
//}

#include "StreamingBcoInfov1.h"

#include <phool/phool.h>

#include <iostream>

void StreamingBcoInfov1::Reset()
{
  set_bco(0);
  set_evtno(0);
  set_usable_bco_tag(false);
  set_bco_streaming_window(std::make_pair(0, 0));
  return;
}

void StreamingBcoInfov1::identify(std::ostream& os) const
{
  os << "identify yourself: I am a StreamingBcoInfov1 Object\n";  return;
}

//int StreamingBcoInfov1::isValid() const
//{
//  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
//  return 0;
//}

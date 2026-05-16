#include "StreamingLumiInfov1.h"

#include <phool/phool.h>

#include <iostream>

void StreamingLumiInfov1::Reset()
{
  // Double check that this is only called once per run!! Or... it should just be called in the InitRun hook?
  set_lumi_raw(0.);
  set_lumi_live(0.);
  set_lumi_scaled(0.);

  return;
}

void StreamingLumiInfov1::identify(std::ostream& os) const
{
  os << "identify yourself: I am a StreamingLumiInfov1 Object\n";  return;
}

//int StreamingLumiInfov1::isValid() const
//{
//  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
//  return 0;
//}

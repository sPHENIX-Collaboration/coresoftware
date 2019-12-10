#include "recoConsts.h"

recoConsts* recoConsts::__instance = nullptr;

void recoConsts::Print() const
{
  // methods from PHFlag
  PrintCharFlags();
  PrintFloatFlags();
  PrintIntFlags();

  return;
}

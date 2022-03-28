#include "recoConsts.h"

recoConsts* recoConsts::__instance = nullptr;

void recoConsts::Print() const
{
  // methods from PHFlag
  PrintStringFlags();
  PrintFloatFlags();
  PrintIntFlags();
  Printuint64Flags();

  return;
}

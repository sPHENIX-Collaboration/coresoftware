#include "recoConsts.h"

using namespace std;

recoConsts* recoConsts::__instance = NULL;

void 
recoConsts::Print() const
{
  // methods from PHFlag
  PrintCharFlags();
  PrintFloatFlags();
  PrintIntFlags();

  return;
}

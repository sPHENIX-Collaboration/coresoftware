#include "VariableArrayUtils.h"

#include <half/half.h>

short VariableArrayUtils::FloatToShortBits(const float rval)
{
  half ftoi(rval);
  return ftoi.bits();
}

float VariableArrayUtils::ShortBitsToFloat(const short ival)
{
  half halfvar;
  halfvar.setBits(ival);
  return halfvar;
}

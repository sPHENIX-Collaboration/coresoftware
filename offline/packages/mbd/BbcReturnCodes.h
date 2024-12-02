// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef BBC_BBCRETURNCODES_H
#define BBC_BBCRETURNCODES_H

#include <limits>

namespace BbcReturnCodes
{
  const short BBC_INVALID_SHORT = std::numeric_limits<short>::min();  //-9999;
  const int BBC_INVALID_INT = std::numeric_limits<int>::min();        //-9999;
  const float BBC_INVALID_FLOAT = std::numeric_limits<float>::quiet_NaN();
}  // namespace BbcReturnCodes

#endif

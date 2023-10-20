
#ifndef MBD_MBDRETURNCODES_H
#define MBD_MBDRETURNCODES_H

#include <limits>

namespace MbdReturnCodes
{
  const short MBD_INVALID_SHORT = std::numeric_limits<short>::min();  //-9999;
  const int MBD_INVALID_INT = std::numeric_limits<int>::min();        //-9999;
  const float MBD_INVALID_FLOAT = std::numeric_limits<float>::quiet_NaN();
}  // namespace MbdReturnCodes

#endif

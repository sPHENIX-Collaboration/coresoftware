// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLUTILS_H
#define FUN4ALL_FUN4ALLUTILS_H

#include <cstdint>
#include <string>
#include <utility>  // for pair

namespace Fun4AllUtils
{
  struct memoryvals
  {
    uint64_t arena{0};
    uint64_t uordblks{0};
    uint64_t fordblks{0};
    uint64_t hblkhd{0};
    uint64_t VmHWM{0};
    uint64_t VmRSS{0};
  };

  std::pair<int, int> GetRunSegment(const std::string &filename);
  void GetMemory(memoryvals &ret);
}  // namespace Fun4AllUtils

#endif

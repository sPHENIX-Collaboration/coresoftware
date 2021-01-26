// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLUTILS_H
#define FUN4ALL_FUN4ALLUTILS_H

#include <string>
#include <utility>  // for pair

namespace Fun4AllUtils
{
  std::pair<int, int> GetRunSegment(const std::string &filename);
}

#endif

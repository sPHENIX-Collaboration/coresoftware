// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHOOL_PHUTILS_H
#define PHOOL_PHUTILS_H

#include <string>

namespace PHUtils
{
  std::string CreateReproducibleTFileName(const std::string &filename);
}

#endif

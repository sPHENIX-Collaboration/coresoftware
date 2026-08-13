#include "PHUtils.h"

#include <phool/phool.h>

#include <filesystem>
#include <iostream>

std::string PHUtils::CreateReproducibleTFileName(const std::string &filename)
{
  std::string outfilename = filename;
  if (filename.empty())
  {
    std::cout << PHWHERE << " called with empty filename string, returning empty string" << std::endl;
    return outfilename;
  }
  std::filesystem::path p = filename;
  outfilename = outfilename + std::string("?reproducible=") + std::string(p.filename());
  return outfilename;
}

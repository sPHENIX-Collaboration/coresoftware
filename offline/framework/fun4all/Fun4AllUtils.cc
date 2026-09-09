#include "Fun4AllUtils.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <malloc.h>
#include <algorithm>  // for max
#include <fstream>
#include <iostream>
#include <vector>

// relying on our standard filenames ...-<runnumber>-<segment>.<ext>
// extract run number and segment number from filename
std::pair<int, int>
Fun4AllUtils::GetRunSegment(const std::string& filename)
{
  int runnumber = 0;
  int segment = -9999;
  boost::char_separator<char> sep("-.");
  boost::tokenizer<boost::char_separator<char> > tok(filename, sep);
  // tokenizer does not have reverse iterator, so fill it in vector
  // and reverse iterate on vector
  std::vector<std::string> tokens;
  for (const auto& t : tok)
  {
    tokens.push_back(t);
  }
  tokens.pop_back();  // remove the file extension
  // try to extract segment number
  try
  {
    segment = boost::lexical_cast<int>((*(tokens.rbegin())));
  }
  catch (boost::bad_lexical_cast const&)
  {
    std::cout << "Cannot extract segment number from filename "
              << filename << std::endl;
    std::cout << "Segment string after parsing: input string "
              << *(tokens.rbegin())
              << " is not valid segment number" << std::endl;
    std::cout << "filename " << filename << " not standard -runnumber-segment.ext"
              << std::endl;
    std::cout << "using " << segment << " as segment number" << std::endl;
  }
  tokens.pop_back();  // remove the segment number
  // try to extract run number
  try
  {
    runnumber = boost::lexical_cast<int>((*(tokens.rbegin())));
  }
  catch (boost::bad_lexical_cast const&)
  {
    std::cout << "Cannot extract run number from filename "
              << filename << std::endl;
    std::cout << "Segment string after parsing: input string "
              << *(tokens.rbegin())
              << " is not valid run number" << std::endl;
    std::cout << "filename " << filename << " not standard -runnumber-segment.ext"
              << std::endl;
    std::cout << "returning " << runnumber << " as run number" << std::endl;
  }
  return std::make_pair(runnumber, segment);
}

void Fun4AllUtils::GetMemory(memoryvals& ret)
{
  const struct mallinfo2 m = mallinfo2();
  ret.arena = m.arena;
  ret.uordblks = m.uordblks;
  ret.fordblks = m.fordblks;
  ret.hblkhd = m.hblkhd;
  std::ifstream status("/proc/self/status");
  std::string line;
  static uint64_t save_VmHWM = 0;
  while (std::getline(status, line))
  {
    if (line.rfind("VmRSS:", 0) == 0)
    {
      std::istringstream fields(line);
      std::string label;
      uint64_t kib = 0;
      std::string unit;

      fields >> label >> kib >> unit;
      ret.VmRSS = kib;
    }
    else if (line.rfind("VmHWM:", 0) == 0)
    {
      std::istringstream fields(line);
      std::string label;
      uint64_t kib = 0;
      std::string unit;

      fields >> label >> kib >> unit;
      ret.VmHWM = kib;
      if (save_VmHWM > kib)
      {
        std::cout << "VmHWM has decreased from " << save_VmHWM << " to " << kib
                  << std::endl;
        std::cout << "line: " << line << std::endl;
      }
      save_VmHWM = kib;
    }
  }
}

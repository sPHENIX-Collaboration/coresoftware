#include "Fun4AllUtils.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm>  // for max
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
  for (auto& t : tok)
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

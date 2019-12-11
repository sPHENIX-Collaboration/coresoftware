#include "Fun4AllUtils.h"

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700)
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <iostream>
#include <vector>

using namespace std;

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
  vector<string> tokens;
  BOOST_FOREACH (string t, tok)
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
    cout << "Cannot extract segment number from filename "
         << filename << endl;
    cout << "Segment string after parsing: input string "
         << *(tokens.rbegin())
         << " is not valid segment number" << endl;
    cout << "filename " << filename << " not standard -runnumber-segment.ext"
         << endl;
    cout << "using " << segment << " as segment number" << endl;
  }
  tokens.pop_back();  // remove the segment number
  // try to extract run number
  try
  {
    runnumber = boost::lexical_cast<int>((*(tokens.rbegin())));
  }
  catch (boost::bad_lexical_cast const&)
  {
    cout << "Cannot extract run number from filename "
         << filename << endl;
    cout << "Segment string after parsing: input string "
         << *(tokens.rbegin())
         << " is not valid run number" << endl;
    cout << "filename " << filename << " not standard -runnumber-segment.ext"
         << endl;
    cout << "returning " << runnumber << " as run number" << endl;
  }
  return make_pair(runnumber, segment);
}

#include "PHSartreGenTrigger.h"

#include <sstream>

using namespace std;

//__________________________________________________________
PHSartreGenTrigger::PHSartreGenTrigger(const std::string &name)
  : m_Verbosity(0)
  , m_Name(name)
{
}

//__________________________________________________________
PHSartreGenTrigger::~PHSartreGenTrigger() {}

std::vector<int> PHSartreGenTrigger::convertToInts(std::string s)
{
  vector<int> theVec;
  stringstream ss(s);
  int i;
  while (ss >> i)
  {
    theVec.push_back(i);
    if (ss.peek() == ',' ||
        ss.peek() == ' ' ||
        ss.peek() == ':' ||
        ss.peek() == ';') ss.ignore();
  }

  return theVec;
}

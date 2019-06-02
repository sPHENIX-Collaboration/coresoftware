#include "PHPy6GenTrigger.h"

using namespace std;

//__________________________________________________________
PHPy6GenTrigger::PHPy6GenTrigger(const std::string &name)
  : _verbosity(0)
  , _name(name)
{
}

//__________________________________________________________
PHPy6GenTrigger::~PHPy6GenTrigger() {}

std::vector<int> PHPy6GenTrigger::convertToInts(std::string s)
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

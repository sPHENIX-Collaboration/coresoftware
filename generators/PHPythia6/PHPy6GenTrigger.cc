#include "PHPy6GenTrigger.h"

#include <sstream>

//__________________________________________________________
PHPy6GenTrigger::PHPy6GenTrigger(const std::string &name)
  : m_Name(name)
{
}

//__________________________________________________________
PHPy6GenTrigger::~PHPy6GenTrigger() {}

std::vector<int> PHPy6GenTrigger::convertToInts(std::string s)
{
  std::vector<int> theVec;
  std::stringstream ss(s);
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

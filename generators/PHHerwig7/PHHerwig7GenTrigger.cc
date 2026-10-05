#include "PHHerwig7GenTrigger.h"

#include <sstream>

//__________________________________________________________
PHHerwig7GenTrigger::PHHerwig7GenTrigger(const std::string& name)
  : m_Verbosity(0)
  , m_Name(name)
{
}

//__________________________________________________________
std::vector<int> PHHerwig7GenTrigger::convertToInts(const std::string& s)
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
        ss.peek() == ';')
    {
      ss.ignore();
    }
  }

  return theVec;
}

#include "Fun4AllBase.h"

#include <iostream>
#include <string>

Fun4AllBase::Fun4AllBase(const std::string& name)
  : m_ThisName(name)
{
  if (name.empty())
  {
    std::cout << "Fun4AllBase::Fun4AllBase: No empty strings as Object Name" << std::endl;
    std::cout << "You likely create a module with an empty name in your macro" << std::endl;
    std::cout << "Since it does not have a name I cannot tell you what it is, but this message is from its ctor, so it happens when you do a new XXX() in your macro" << std::endl;
    exit(1);
  }
  return;
}

Fun4AllBase::~Fun4AllBase()
{
  if (m_Verbosity >= VERBOSITY_MORE)
  {
    std::cout << "Deleting " << m_ThisName << std::endl;
  }
  return;
}

void Fun4AllBase::Print(const std::string& /*what*/) const
{
  std::cout << Name() << " did not implement Print method" << std::endl;
  return;
}

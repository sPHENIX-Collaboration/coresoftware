#include "Fun4AllBase.h"

#include <iostream>
#include <string>

Fun4AllBase::Fun4AllBase(const std::string& name)
  : m_ThisName(name)
{
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

#include "Fun4AllBase.h"

#include <iostream>

using namespace std;

Fun4AllBase::Fun4AllBase(const string &name):
  m_ThisName(name),
  m_Verbosity(VERBOSITY_QUIET),
  ThisName(name),
  verbosity(VERBOSITY_QUIET)  
{
  return;
}

Fun4AllBase::~Fun4AllBase()
{
  if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "Deleting " << Name () << endl;
    }
  return;
}

void Fun4AllBase::Print(const string& /*what*/) const
{
  cout << Name() << " did not implement Print method" << endl;
  return;
}

void
Fun4AllBase::Name(const std::string &name)
{
  m_ThisName = name;
  ThisName = name;
}

const string
Fun4AllBase::Name() const
{
  m_ThisName = ThisName;
  return m_ThisName;
}

void
Fun4AllBase::Verbosity(const int ival)
{
  m_Verbosity = ival;
  verbosity = ival;
  return;
}

void
Fun4AllBase::Verbosity(enu_Verbosity ival) 
{
  m_Verbosity = ival;
  verbosity = ival;
  return;
}

int
Fun4AllBase::Verbosity() const
{
  m_Verbosity = verbosity;
  return m_Verbosity;
}


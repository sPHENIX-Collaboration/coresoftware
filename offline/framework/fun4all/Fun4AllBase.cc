#include "Fun4AllBase.h"

#include <iostream>

using namespace std;

Fun4AllBase::Fun4AllBase(const string& name)
  : m_ThisName(name)
  , m_Verbosity(VERBOSITY_QUIET)
{
  return;
}

Fun4AllBase::~Fun4AllBase()
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "Deleting " << Name() << endl;
  }
  return;
}

void Fun4AllBase::Print(const string& /*what*/) const
{
  cout << Name() << " did not implement Print method" << endl;
  return;
}

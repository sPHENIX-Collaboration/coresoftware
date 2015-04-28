#include "Fun4AllBase.h"

#include <iostream>

using namespace std;

Fun4AllBase::Fun4AllBase(const string &name):
  ThisName(name),
  verbosity(0)
{
  return;
}

Fun4AllBase::~Fun4AllBase()
{
  if (verbosity > 0)
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


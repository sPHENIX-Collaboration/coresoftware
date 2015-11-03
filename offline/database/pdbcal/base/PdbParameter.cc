#include "PdbParameter.hh"
#include <phool/phool.h>

#include <iostream>

using namespace std;

PdbParameter::PdbParameter()
{
  thePar = 0;
  setName("UNKNOWN");
}

PdbParameter::PdbParameter(const float value)
{
  thePar = value;
  setName("UNKNOWN");
}

PdbParameter::PdbParameter(const float value, const char *name)
{
  thePar = value;
  setName(name);
}

void PdbParameter::print() const
{
  cout << theName << ": " << thePar << endl;
}


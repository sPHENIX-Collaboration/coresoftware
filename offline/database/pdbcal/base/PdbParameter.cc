#include "PdbParameter.h"

#include <cmath>
#include <iostream>

using namespace std;

PdbParameter::PdbParameter()
  : thePar(NAN)
{
}

PdbParameter::PdbParameter(const double value, const string &name)
  : thePar(value)
  , theName(name)
{
}

void PdbParameter::print() const
{
  cout << theName << ": " << thePar << endl;
}

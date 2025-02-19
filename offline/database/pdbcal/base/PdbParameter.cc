#include "PdbParameter.h"

#include <iostream>

PdbParameter::PdbParameter(const double value, const std::string &name)
  : thePar(value)
  , theName(name)
{
}

void PdbParameter::print() const
{
  std::cout << theName << ": " << thePar << std::endl;
}

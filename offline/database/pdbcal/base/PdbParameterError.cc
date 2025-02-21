//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 2004
//
//  Implementation of class PdbParameter
//
//  Author: Cesar

#include "PdbParameterError.h"

#include <cmath>
#include <iostream>

PdbParameterError::PdbParameterError(const double value, const double error, const std::string &name)
  : PdbParameter(value, name)
  , theParError(error)
{
}

void PdbParameterError::print() const
{
  std::cout << getName() << ": " << getParameter() << " +/- " << getParameterError() << std::endl;
}

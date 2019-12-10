//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 2004
//
//  Implementation of class PdbParameter
//
//  Author: Cesar

#include "PdbParameterError.h"

#include <cmath>
#include <iostream>

using namespace std;

PdbParameterError::PdbParameterError()
  : theParError(NAN)
{
}

PdbParameterError::PdbParameterError(const double value, const double error, const string &name)
  : PdbParameter(value, name)
  , theParError(error)
{
}

void PdbParameterError::print() const
{
  cout << theName << ": " << thePar << " +/- " << theParError << endl;
}

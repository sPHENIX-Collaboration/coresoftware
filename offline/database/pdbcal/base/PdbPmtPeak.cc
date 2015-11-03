//-----------------------------------------------------------------------------
//  $header$
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbPmtPeak
//
//  Author: ohnishi
//-----------------------------------------------------------------------------
#include "PdbPmtPeak.hh"

#include <iostream>

PdbPmtPeak::PdbPmtPeak()
{
  PeakChannel = -9999.9;
  Deviation = -9999.9;
  Status = -9999;
}

void
PdbPmtPeak::print() const
{
  std::cout << std::endl;
  std::cout << "Peak Channel : " << PeakChannel << std::endl;
  std::cout << "Deviation    : " << Deviation << std::endl;
  std::cout << "Status       : " << Status << std::endl;
}

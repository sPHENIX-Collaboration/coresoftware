//  Implementation of class PdbPmtFitPar
//
//  Author: ohnishi

#include <iostream>
#include "PdbPmtFitPar.hh"

using namespace std;

PdbPmtFitPar::PdbPmtFitPar()
{
  Par0 = 0.0;
  Par1 = 0.0;
  Par2 = 0.0;
  Par3 = 0.0;
  Par4 = 0.0;
  Chi2 = 0.0;
  Status = 0;
}

PdbPmtFitPar::~PdbPmtFitPar()
{}

PdbPmtFitPar& PdbPmtFitPar::operator=(const  PdbPmtFitPar&p)
{
  Par0 = p.Par0;
  Par1 = p.Par1;
  Par2 = p.Par2;
  Par3 = p.Par3;
  Par4 = p.Par4;
  Chi2 = p.Chi2;
  Status = p.Status;
  return *this;
}

void
PdbPmtFitPar::print() const
{
  cout << "nPmtFitPar = "
       << Par0 << " " << Par1 << " " 
       << Par2 << " " << Par3 << " " 
       << Par4 
       << " Chi2 = " << Chi2 << endl;
}

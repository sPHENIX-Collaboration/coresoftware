//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbBankID2.cc,v 1.1 2010/09/10 17:00:43 irina Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbBankID
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PdbBankID2.h"
#include "PHString.h"
#include <iostream>

ClassImp(PdbBankID2);

PdbBankID2::PdbBankID2(const int val)
{
   bankID = val;
}

void PdbBankID2::print() const
{
  std::cout << "BankID of PdbBankID2= " << bankID  << std::endl;
}

int operator == (const PdbBankID2 & a, const PdbBankID2 & b)
{
   return (a.bankID == b.bankID);
}

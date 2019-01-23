//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbBankID.cc,v 1.4 2006/01/10 06:39:47 pinkenbu Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbBankID
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PdbBankID.h"
#include <iostream>

PdbBankID::PdbBankID(const int val)
{
  bankID = val;
}

void PdbBankID::print() const
{
  std::cout << "BankID = " << bankID << std::endl;
}

int operator==(const PdbBankID& a, const PdbBankID& b)
{
  return (a.bankID == b.bankID);
}

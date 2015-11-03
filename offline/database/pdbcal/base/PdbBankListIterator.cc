//-----------------------------------------------------------------------------
//  $header$
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbBankListIterator
//
//  Author: messer
//-----------------------------------------------------------------------------
#include "PdbBankListIterator.hh"

#include <iostream>

PdbBankListIterator::PdbBankListIterator()
{
}

PdbBankListIterator::PdbBankListIterator(PdbBankList & bankList) : PHPointerListIterator<PdbCalBank>(bankList)
{
}

PdbBankListIterator::~PdbBankListIterator()
{
}

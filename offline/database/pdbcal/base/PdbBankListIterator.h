//-----------------------------------------------------------------------------
//
//  The pdbcal package
//
//  Declaration of class PdbBankListIterator
//
//  Purpose: iterator class for PdbBankList
//
//  Description:
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef PDBCAL_BASE_PDBBANKLISTITERATOR_H
#define PDBCAL_BASE_PDBBANKLISTITERATOR_H

#include "PdbBankList.h"
#include "PdbCalBank.h"

#include <phool/PHPointerListIterator.h>

class PdbBankListIterator : public PHPointerListIterator<PdbCalBank>
{
 public:
  PdbBankListIterator(PdbBankList &bankList)
    : PHPointerListIterator<PdbCalBank>(bankList)
  {
  }

  ~PdbBankListIterator() override {}
 private: 
  PdbBankListIterator() = delete;
};

#endif  // PDBCAL_BASE_PDBBANKLISTITERATOR_H

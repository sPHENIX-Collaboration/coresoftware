//-----------------------------------------------------------------------------
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbBankListIterator
//
//  Purpose: iterator class for PdbBankList
//
//  Description:
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef __PDBBANKLISTITERATOR_HH__
#define __PDBBANKLISTITERATOR_HH__

#include "PdbBankList.h"
#include "PdbCalBank.h"

#include <phool/PHPointerListIterator.h>

class PdbBankListIterator : public PHPointerListIterator<PdbCalBank> {
public:
   PdbBankListIterator(PdbBankList &);
   ~PdbBankListIterator();
   
protected:
   PdbBankListIterator();

};

#endif /* __PDBBANKLISTITERATOR_HH__ */

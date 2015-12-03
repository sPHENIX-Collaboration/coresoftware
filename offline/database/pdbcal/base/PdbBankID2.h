//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbBankID2.hh,v 1.1 2010/09/10 16:59:50 irina Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbBankID
//
//  Purpose: id number for a bank, derived from a string
//
//  Description: The string should follow the ONCS naming convention
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef __PDBBANKID2_HH__
#define __PDBBANKID2_HH__

#include <TObject.h>

class PdbBankID2 : public TObject {
public:
   PdbBankID2(const int val = 0);
   virtual ~PdbBankID2(){}

   void print() const;

   int  getInternalValue() const { return bankID; }
   void setInternalValue(const int val) { bankID = val; }
   
   friend int operator == (const PdbBankID2 &, const PdbBankID2 &);
   
private:
  int bankID;

  ClassDef(PdbBankID2, 1);

};

#endif /* __PDBBANKID2_HH__ */

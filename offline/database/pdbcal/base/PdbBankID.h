//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbBankID.hh,v 1.6 2006/02/04 16:04:06 pinkenbu Exp $
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
#ifndef __PDBBANKID_HH__
#define __PDBBANKID_HH__

class PdbBankID {
public:
   PdbBankID(const int val = 0);
   virtual ~PdbBankID(){}

   void print() const;

   int  getInternalValue() const { return bankID; }
   void setInternalValue(const int val) { bankID = val; }
   
   friend int operator == (const PdbBankID &, const PdbBankID &);
   
private:
  int bankID;

};

#endif /* __PDBBANKID_HH__ */

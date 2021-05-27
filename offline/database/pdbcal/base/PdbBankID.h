//-----------------------------------------------------------------------------
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
#ifndef PDBCAL_BASE_PDBBANKID_H
#define PDBCAL_BASE_PDBBANKID_H

#include <phool/PHObject.h>

class PdbBankID : public PHObject {
public:
   PdbBankID(const int val = 0);
   ~PdbBankID() override{}

   void print() const;

   int  getInternalValue() const { return bankID; }
   void setInternalValue(const int val) { bankID = val; }
   
   friend int operator == (const PdbBankID &, const PdbBankID &);
   
private:
  int bankID;

  ClassDefOverride(PdbBankID, 1)

};

#endif /* PDBCAL_BASE_PDBBANKID_H */

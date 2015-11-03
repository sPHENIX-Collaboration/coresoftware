//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbCalHeader.hh,v 1.2 2008/05/31 21:40:23 pinkenbu Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbCalHeader
//
//  Purpose: Stores the validity range and ID of a bank
//
//  Description:
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef __PDBCALHEADER_HH__
#define __PDBCALHEADER_HH__

#include "Pdb.hh"
#include "PdbBankID.hh"
#include <phool/PHTimeStamp.h>

#include <string>

class PdbCalHeader {
public:
   PdbCalHeader();
   PdbCalHeader(PdbBankID, const char *, PHTimeStamp &, PHTimeStamp &);
   PdbCalHeader( const PdbCalHeader& );

  virtual ~PdbCalHeader( ){}

   const PdbBankID&   getBankID()       const { return bankID; }
   const PHTimeStamp& getInsertTime()   const { return insertTime; }
   const PHTimeStamp& getStartValTime() const { return startValTime; }
   const PHTimeStamp& getEndValTime()   const { return endValTime; }
   const std::string     getDescription()  const { return description; }
   const std::string     getUserName()     const { return userName; }
   
   void setBankID(const PdbBankID & val)         { bankID = val; }
   void setInsertTime(const PHTimeStamp & val)   { insertTime = val; }
   void setStartValTime(const PHTimeStamp & val) { startValTime = val; }
   void setEndValTime(const PHTimeStamp & val)   { endValTime = val; }
  void setDescription(const std::string &val) {description = val;}
  void setUserName(const std::string &val) {userName = val;}
   
   void print() const;

private:
  PdbBankID   bankID;
  PHTimeStamp insertTime;
  PHTimeStamp startValTime;
  PHTimeStamp endValTime;
  
  std::string description;	// Description of calibration.
  std::string userName;		// The user who made the bank.

};

#endif /* __PDBCALHEADER_HH__ */

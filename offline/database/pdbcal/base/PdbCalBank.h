 //  Declaration of class PdbCalBank
//  Purpose: Calibration bank base class
//  Author: Matthias Messer

#ifndef __PDBCALBANK_HH__
#define __PDBCALBANK_HH__

#include "PdbBankID.h"
#include "PdbBankID2.h"
#include "phool/PHTimeStamp.h"
#include "PHString.h"
#ifndef __CINT__
#include <cstddef>
#endif

#include <TObject.h>

class PdbCalChan;
class PHTimeStamp;
class PdbBankID;
class PdbBankID2;

class PdbCalBank : public  TObject 
{
public:
   PdbCalBank();
   virtual ~PdbCalBank();
   virtual PdbCalBank* clone() const = 0;

   virtual void printHeader() const = 0;
   virtual void print() = 0;
   virtual void printEntry(size_t) = 0;
   
   virtual size_t         getLength() = 0;
   virtual PdbCalChan &   getEntry(size_t) = 0;
   virtual void setLength(size_t val) = 0;

   //
   // Access functions for the header
   //
   virtual PdbBankID   getBankID()       const = 0;
   virtual PdbBankID2  getBankID2()      const = 0;
   virtual PHTimeStamp getInsertTime()   const = 0;
   virtual PHTimeStamp getStartValTime() const = 0;
   virtual PHTimeStamp getEndValTime()   const = 0;
   virtual PHString    getDescription()  const = 0;
   virtual PHString    getUserName()     const = 0;
   
   virtual void setBankID(const PdbBankID &)         = 0; 
   virtual void setBankID2(const PdbBankID2 &)       = 0;
   virtual void setInsertTime(const PHTimeStamp &)   = 0;
   virtual void setStartValTime(const PHTimeStamp &) = 0;
   virtual void setEndValTime(const PHTimeStamp &)   = 0;
   virtual void setDescription(const PHString &)     = 0;
   virtual void setUserName(const PHString &)        = 0;
   
   virtual int isValid(const PHTimeStamp &) const = 0;

  ClassDef(PdbCalBank,2);
};

#endif /* __PDBCALBANK_HH__ */

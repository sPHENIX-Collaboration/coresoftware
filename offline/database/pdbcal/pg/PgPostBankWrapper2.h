#ifndef __PGPOSTBANKWRAPPER2_HH__
#define __PGPOSTBANKWRAPPER2_HH__

#include "PgPostCalBank.h"

#include <phool/PHTimeStamp.h>

#include <iostream>

class PgPostBankWrapper2 : public PgPostCalBank {
public:
  PgPostBankWrapper2 ();
  PgPostBankWrapper2 (PdbCalBank *b);
  virtual ~PgPostBankWrapper2 ();
  virtual PgPostCalBank * clone() const { return new PgPostBankWrapper2(*this); }
  
  void printHeader () const;
  void print (){ bank->print(); }
  void printEntry (size_t s) { bank->printEntry(s); }
   
  size_t getLength () {return bank->getLength(); }
  PdbCalChan& getEntry (size_t pos) { return bank->getEntry(pos); } 
  void setLength (size_t len) { bank->setLength(len); }
  virtual bool commit();
  virtual bool commit_rid(int rid,long it,long st,long et);

  PdbBankID2    getBankID2()    const { return bankID2; }
  PHTimeStamp getInsertTime()   const { return insertTime; }
  PHTimeStamp getStartValTime() const { return startValTime; }
  PHTimeStamp getEndValTime()   const { return endValTime; }
  std::string    getDescription()  const { return description; }
  std::string    getUserName()     const { return userName; }
  std::string    getTableName()    const { return tableName; }
 
  void setBankID2(const PdbBankID2 & val)          { bankID2 = val; }
  void setInsertTime(const PHTimeStamp & val)   { insertTime = val; }
  void setStartValTime(const PHTimeStamp & val) { startValTime = val; }
  void setEndValTime(const PHTimeStamp & val)   { endValTime = val; }
  void setDescription(const std::string & val)     { strcpy(description, val.c_str()); }  
  void setUserName(const std::string & val)        { strcpy(userName, val.c_str()); }
  void setTableName(const std::string & val)       { strcpy(tableName, val.c_str()); }

  PdbCalBank * getBank() { return  bank; }
  virtual int isValid (const PHTimeStamp &) const { return 0; }

private:

  PdbBankID2     bankID2;
  PHTimeStamp  insertTime;
  PHTimeStamp  startValTime;
  PHTimeStamp  endValTime;
  char         description[240];
  char         userName[200];
  char         tableName[400];

  PdbCalBank * bank;

  ClassDef(PgPostBankWrapper2, 1);

};

#endif /* __PGPOSTBANKWRAPPER2_HH__ */

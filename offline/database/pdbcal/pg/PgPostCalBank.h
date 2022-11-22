// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTCALBANK_H
#define PDBCALPG_PGPOSTCALBANK_H

#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbCalBank.h>

#include <phool/PHTimeStamp.h>

#include <TObject.h>

#include <cstring>
#include <iostream>

class PgPostCalBank : public PdbCalBank
{
 public:
  PgPostCalBank() {}
  virtual ~PgPostCalBank() {}

  virtual void printHeader() const { std::cout << "I'm PgPostCalBank" << std::endl; }
  virtual void printEntry(size_t) = 0;
  virtual void print() = 0;

  //  virtual bool commit() = 0;
  virtual size_t getLength() = 0;
  virtual PdbCalChan& getEntry(size_t) = 0;
  virtual void setLength(size_t val) = 0;

  virtual PdbBankID getBankID() const { return 0; }
  virtual PHTimeStamp getInsertTime() const { return PHTimeStamp((time_t) 0); }
  virtual PHTimeStamp getStartValTime() const { return PHTimeStamp((time_t) 0); }
  virtual PHTimeStamp getEndValTime() const { return PHTimeStamp((time_t) 0); }
  virtual std::string getDescription() const { return ""; }
  virtual std::string getUserName() const { return ""; }
  virtual std::string getTableName() const { return ""; }

  virtual void setBankID(const PdbBankID& /*val*/) {}
  virtual void setInsertTime(const PHTimeStamp& /*val*/) {}
  virtual void setStartValTime(const PHTimeStamp& /*val*/) {}
  virtual void setEndValTime(const PHTimeStamp& /*val*/) {}
  virtual void setDescription(const std::string& /*val*/) {}
  virtual void setUserName(const std::string& /*val*/) {}
  virtual void setTableName(const std::string& /*val*/) {}

  virtual int isValid(const PHTimeStamp&) const { return 0; }

  ClassDefOverride(PgPostCalBank, 2);
};

#endif

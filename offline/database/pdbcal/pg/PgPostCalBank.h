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

  void printHeader() const override { std::cout << "I'm PgPostCalBank" << std::endl; }
  void printEntry(size_t) override = 0;
  void print() override = 0;

  //  virtual bool commit() = 0;
  size_t getLength() override = 0;
  PdbCalChan& getEntry(size_t) override = 0;
  void setLength(size_t val) override = 0;

  PdbBankID getBankID() const override { return 0; }
  PHTimeStamp getInsertTime() const override { return PHTimeStamp((time_t) 0); }
  PHTimeStamp getStartValTime() const override { return PHTimeStamp((time_t) 0); }
  PHTimeStamp getEndValTime() const override { return PHTimeStamp((time_t) 0); }
  std::string getDescription() const override { return ""; }
  std::string getUserName() const override { return ""; }
  virtual std::string getTableName() const { return ""; }

  void setBankID(const PdbBankID& /*val*/) override {}
  void setInsertTime(const PHTimeStamp& /*val*/) override {}
  void setStartValTime(const PHTimeStamp& /*val*/) override {}
  void setEndValTime(const PHTimeStamp& /*val*/) override {}
  void setDescription(const std::string& /*val*/) override {}
  void setUserName(const std::string& /*val*/) override {}
  virtual void setTableName(const std::string& /*val*/) {}

  int isValid(const PHTimeStamp&) const override { return 0; }

  ClassDefOverride(PgPostCalBank, 2);
};

#endif

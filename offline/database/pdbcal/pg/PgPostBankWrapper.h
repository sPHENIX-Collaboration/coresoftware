// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTBANKWRAPPER_H
#define PDBCALPG_PGPOSTBANKWRAPPER_H

#include "PgPostCalBank.h"

#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbCalBank.h>

#include <phool/PHTimeStamp.h>

#include <cstddef>
#include <string>

class PdbCalChan;
class PHObject;

class PgPostBankWrapper : public PgPostCalBank
{
 public:
  PgPostBankWrapper();
  PgPostBankWrapper(PdbCalBank *b);
  ~PgPostBankWrapper() override;
  PHObject *CloneMe() const override { return new PgPostBankWrapper(*this); }

  void printHeader() const override;
  void print() override { bank->print(); }
  void printEntry(size_t s) override { bank->printEntry(s); }

  size_t getLength() override { return bank->getLength(); }
  PdbCalChan &getEntry(size_t pos) override { return bank->getEntry(pos); }
  void setLength(size_t len) override { bank->setLength(len); }
  virtual bool commit();

  PdbBankID getBankID() const override { return bankID; }
  PHTimeStamp getInsertTime() const override { return insertTime; }
  PHTimeStamp getStartValTime() const override { return startValTime; }
  PHTimeStamp getEndValTime() const override { return endValTime; }
  std::string getDescription() const override { return description; }
  std::string getUserName() const override { return userName; }
  std::string getTableName() const override { return tableName; }

  void setBankID(const PdbBankID &val) override { bankID = val; }
  void setInsertTime(const PHTimeStamp &val) override { insertTime = val; }
  void setStartValTime(const PHTimeStamp &val) override { startValTime = val; }
  void setEndValTime(const PHTimeStamp &val) override { endValTime = val; }
  void setDescription(const std::string &val) override { description = val; }
  void setUserName(const std::string &val) override { userName = val; }
  void setTableName(const std::string &val) override { tableName = val; }

  PdbCalBank *getBank() { return bank; }
  int isValid(const PHTimeStamp &) const override { return 0; }

 private:
  PdbBankID bankID;
  PHTimeStamp insertTime;
  PHTimeStamp startValTime;
  PHTimeStamp endValTime;
  std::string description;
  std::string userName;
  std::string tableName;

  PdbCalBank *bank;

  ClassDefOverride(PgPostBankWrapper, 1);
};

#endif /* PDBCAL_PG_PGPOSTBANKWRAPPER_H */

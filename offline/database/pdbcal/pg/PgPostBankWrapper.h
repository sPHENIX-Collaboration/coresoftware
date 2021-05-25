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
  virtual ~PgPostBankWrapper();
  virtual PHObject *CloneMe() const { return new PgPostBankWrapper(*this); }

  void printHeader() const;
  void print() { bank->print(); }
  void printEntry(size_t s) { bank->printEntry(s); }

  size_t getLength() { return bank->getLength(); }
  PdbCalChan &getEntry(size_t pos) { return bank->getEntry(pos); }
  void setLength(size_t len) { bank->setLength(len); }
  virtual bool commit();

  PdbBankID getBankID() const { return bankID; }
  PHTimeStamp getInsertTime() const { return insertTime; }
  PHTimeStamp getStartValTime() const { return startValTime; }
  PHTimeStamp getEndValTime() const { return endValTime; }
  std::string getDescription() const { return description; }
  std::string getUserName() const { return userName; }
  std::string getTableName() const { return tableName; }

  void setBankID(const PdbBankID &val) { bankID = val; }
  void setInsertTime(const PHTimeStamp &val) { insertTime = val; }
  void setStartValTime(const PHTimeStamp &val) { startValTime = val; }
  void setEndValTime(const PHTimeStamp &val) { endValTime = val; }
  void setDescription(const std::string &val) { description = val; }
  void setUserName(const std::string &val) { userName = val; }
  void setTableName(const std::string &val) { tableName = val; }

  PdbCalBank *getBank() { return bank; }
  virtual int isValid(const PHTimeStamp &) const { return 0; }

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

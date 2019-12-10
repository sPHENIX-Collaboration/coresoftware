// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTCALBANKITERATOR_H
#define PDBCALPG_PGPOSTCALBANKITERATOR_H

#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbCalBankIterator.h>

#include <phool/PHTimeStamp.h>

#include <ctime>
#include <map>
#include <iostream>
#include <string>

class PdbCalBank;
class PgPostApplication;
class PgPostBankManager;
class TSQLResultSet;
class TSQLStatement;

class PgPostCalBankIterator : public PdbCalBankIterator
{
 public:
  PgPostCalBankIterator(PgPostBankManager& bm);
  virtual ~PgPostCalBankIterator();

  virtual bool init(const std::string& fulldbname, const PdbBankID& bankid);

  virtual bool isValid() const;

  virtual void print(std::ostream& os = std::cout) const;

  virtual void setBankID(const PdbBankID& id);

  virtual void setEndValTimeLimits(const PHTimeStamp& min = PHTimeStamp(0),
                                   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));

  virtual void setInsertTimeLimits(const PHTimeStamp& min = PHTimeStamp(0),
                                   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));

  virtual void setStartValTimeLimits(const PHTimeStamp& min = PHTimeStamp(0),
                                     const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));

  virtual PdbCalBank* next();

 private:
  class ValPeriod
  {
   public:
    ValPeriod(const PHTimeStamp& min = PHTimeStamp(0),
              const PHTimeStamp& max = PHTimeStamp::PHFarFuture)
      : start_(min)
      , end_(max)
    {
    }

    time_t start() const { return start_.getTics(); }
    time_t end() const { return end_.getTics(); }

   private:
    PHTimeStamp start_;
    PHTimeStamp end_;
  };

  typedef std::map<std::string, ValPeriod> TimeMap;
  TimeMap fTimeMap;

  PgPostBankManager& fBM;
  PgPostApplication* fApplication;
  std::string fDBName;
  std::string fTableName;
  bool fIsValid;
  PdbBankID fBankID;
  TSQLStatement* fSQLStatement;
  TSQLResultSet* fSQLResultSet;
};

#endif

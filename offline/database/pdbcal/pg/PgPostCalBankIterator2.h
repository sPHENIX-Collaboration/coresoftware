#ifndef __PGPOSTCALBANKITERATOR2_H_HH
#define __PGPOSTCALBANKITERATOR2_H_HH

#include <pdbcalbase/PdbCalBankIterator.h>
#include <pdbcalbase/PdbBankID2.h>

#include <map>
#include <string>

class PgPostBankManager;
class PgPostApplication;
class TSQLResultSet;
class TSQLStatement;

class PgPostCalBankIterator2 : public PdbCalBankIterator
{
 public:
  
  PgPostCalBankIterator2(PgPostBankManager& bm);
  virtual ~PgPostCalBankIterator2();

  virtual bool init(const char* fulldbname, const PdbBankID2& bankid2);
  
  virtual bool isValid() const;

  virtual void print(std::ostream& os = std::cout) const;

  virtual void setBankID2(const PdbBankID2& id);

  virtual void setEndValTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));
  
  virtual void setInsertTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));
  
  virtual void setStartValTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture));

  virtual PdbCalBank* next();

 private:

  class ValPeriod
  {
  public:

    ValPeriod(const PHTimeStamp& min=PHTimeStamp(0), 
	      const PHTimeStamp& max=PHTimeStamp::PHFarFuture)
      : start_(min), end_(max) {}

    time_t start() const { return start_.getTics(); }
    time_t end() const { return end_.getTics(); }

  private:
    PHTimeStamp start_;
    PHTimeStamp end_;
  };

  typedef std::map<std::string,ValPeriod> TimeMap;
  TimeMap fTimeMap;

  PgPostBankManager& fBM;
  PgPostApplication* fApplication;
  std::string fDBName;
  std::string fTableName;
  bool fIsValid;
  PdbBankID2 fBankID2;
  TSQLStatement* fSQLStatement;
  TSQLResultSet* fSQLResultSet;
};

#endif

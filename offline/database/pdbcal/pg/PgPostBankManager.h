// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTBANKMANAGER_H
#define PDBCALPG_PGPOSTBANKMANAGER_H

#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankManager.h>

#include <phool/PHTimeStamp.h>

#include <map>
#include <set>
#include <string>
#include <ctime> 

class PdbApplication;
class PdbCalBank;
class PdbCalBankIterator;

class PgPostBankManager : public PdbBankManager
{
 public:
  static PgPostBankManager *instance();
  static int Register();
  virtual ~PgPostBankManager();

 protected:
  PgPostBankManager();

 public:
  virtual PdbCalBankIterator *getIterator();
  virtual PdbCalBank *createBank(const std::string &,
                                 PdbBankID,
                                 const std::string &, PHTimeStamp &, PHTimeStamp &, const std::string &);

  // create bank with run number as key
  PdbCalBank *createBank(const int, const std::string &, PdbBankID, const std::string &, const std::string &, const time_t duration = 60);

  // create bank for a given range of run numbers rather than timestamps
  PdbCalBank *createBank(const int, const int, const std::string &, PdbBankID, const std::string &, const std::string &);

  // fetch banks with run number as key
  PdbCalBank *fetchBank(const std::string &, PdbBankID, const std::string &, const int);

  PdbCalBank *fetchClosestBank(const std::string &, PdbBankID, const std::string &, const int);

  //   void fetchAllBanks(PdbBankList &, const std::string &, PdbBankID, const std::string &, const int);
  //   void fetchAllBanks(PdbBankList &, const std::string &, const std::string &, const int);

  // fetch banks with timestamp as key
  PdbCalBank *fetchBank(const std::string &, PdbBankID, const std::string &, const PHTimeStamp &);

  PdbCalBank *fetchClosestBank(const std::string &, PdbBankID, const std::string &, PHTimeStamp &);

  // void fetchAllBanks(PdbBankList &, const std::string &, PdbBankID, const std::string &, PHTimeStamp &);

  // void fetchAllBanks(PdbBankList &, const std::string &, const std::string &, PHTimeStamp &);

  PdbApplication *getApplication();

  void fillCalibObject(PdbCalBank *, const std::string &, PHTimeStamp &) {}

  void GetUsedBankRids(std::map<std::string, std::set<int> > &usedbanks) const;
  void ClearUsedBankRids() { BankRid.clear(); }
  void SetMaxInsertTime(const PHTimeStamp &tMax);

 private:
  static PgPostBankManager *mySpecificCopy;
  std::map<std::string, std::set<int> > BankRid;
  std::string getRealName(const std::string &);
  PHTimeStamp tMaxInsertTime;
};

#endif /* PDBCAL_PG_PGPOSTBANKMANAGER_H */

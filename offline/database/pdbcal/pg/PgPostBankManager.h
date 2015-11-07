// $Id: PgPostBankManager.h,v 1.13 2014/03/21 20:54:44 phnxbld Exp $

#ifndef __PGPOSTBANKMANAGER_HH__
#define __PGPOSTBANKMANAGER_HH__

#include <pdbcalbase/PdbBankManager.h>

#include <map>
#include <set>
#include <string>

class PgPostBankManager : public PdbBankManager {

public:
  static PgPostBankManager *instance();
  static int Register();
  virtual ~PgPostBankManager( );

protected:  
  PgPostBankManager();

public:
 

  virtual PdbCalBankIterator* getIterator();
  virtual PdbCalBank* createBank(const char *,
				 PdbBankID,
				 const char *, PHTimeStamp &,PHTimeStamp &,const char *);


  // create bank with run number as key
  PdbCalBank* createBank(const int, const char *, PdbBankID, const char *, const char *, const time_t duration=60);

  // create bank for a given range of run numbers rather than timestamps
  PdbCalBank* createBank(const int, const int, const char *, PdbBankID, const char *, const char *);

  // fetch banks with run number as key
  PdbCalBank* fetchBank(const char *, PdbBankID, const char *, const int);

  PdbCalBank* fetchClosestBank(const char *, PdbBankID, const char *, const int);


//   void fetchAllBanks(PdbBankList &, const char *, PdbBankID, const char *, const int);
// void fetchAllBanks(PdbBankList &, const char *, PdbBankID2, const char *, const int);

//   void fetchAllBanks(PdbBankList &, const char *, const char *, const int);

  // fetch banks with timestamp as key
  PdbCalBank* fetchBank(const char *, PdbBankID, const char *, const PHTimeStamp &);


  PdbCalBank* fetchClosestBank(const char *, PdbBankID, const char *, PHTimeStamp &);

  // void fetchAllBanks(PdbBankList &, const char *, PdbBankID, const char *, PHTimeStamp &);
  // void fetchAllBanks(PdbBankList &, const char *, PdbBankID2, const char *, PHTimeStamp &);

  // void fetchAllBanks(PdbBankList &, const char *, const char *, PHTimeStamp &);
  
  PdbApplication* getApplication(PHBoolean pJob = False);
  
  void fillCalibObject(PdbCalBank*, const char *,  PHTimeStamp &) {}

  void GetUsedBankRids(std::map<std::string,std::set<int> > &usedbanks) const;
  void ClearUsedBankRids() {BankRid.clear();}
  void SetMaxInsertTime(const PHTimeStamp &tMax);

private:
  static PgPostBankManager *mySpecificCopy;
  std::map<std::string,std::set<int> > BankRid;
  std::string getRealName(const std::string &);
  PHTimeStamp tMaxInsertTime;
};

#endif /* __PGPOSTBANKMANAGER_HH__ */

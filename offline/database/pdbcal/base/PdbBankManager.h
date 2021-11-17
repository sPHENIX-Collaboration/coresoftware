//  Declaration of class PdbBankManager
//  Purpose: Abstract factory class to handle banks
//  Author: Matthias Messer

#ifndef PDBCAL_BASE_PDBBANKMANAGER_H
#define PDBCAL_BASE_PDBBANKMANAGER_H

#include "PdbBankID.h"

#include <phool/PHTimeStamp.h>

#include <ctime>
#include <map>
#include <set>
#include <string>

class PdbCalBank;
class PdbApplication;
class PdbCalBankIterator;

class PdbBankManager 
{

protected:

  PdbBankManager(){}
  virtual ~PdbBankManager();

public:

  static  PdbBankManager *instance();

  /// Get an iterator to loop over banks.
  virtual PdbCalBankIterator* getIterator() = 0;

  virtual PdbCalBank* createBank(const std::string &,
				 PdbBankID,
				 const std::string &,
				 PHTimeStamp &,
				 PHTimeStamp &,
				 const std::string &) = 0;

  // create bank with run number as key
  virtual PdbCalBank* createBank(const int,
				 const std::string &,
				 PdbBankID,
				 const std::string &,
				 const std::string &,
				 const time_t duration=60) = 0;

  // create bank for a given range of run numbers rather than timestamps
  virtual PdbCalBank* createBank(const int,
				 const int,
				 const std::string &,
				 PdbBankID,
				 const std::string &,
				 const std::string &) = 0;
  
  virtual PdbCalBank* fetchBank(const std::string &,
				PdbBankID,
				const std::string &,
				const int) = 0;

  virtual PdbCalBank* fetchClosestBank(const std::string &,
				       PdbBankID,
				       const std::string &,
				       const int) = 0; 
  // virtual void fetchAllBanks(PdbBankList &,
  // 			     const std::string &,
  // 			     PdbBankID,
  // 			     const std::string &,
  // 			     const int) = 0;

  // virtual void fetchAllBanks(PdbBankList &,
  // 			     const std::string &,
  // 			     const std::string &,
  // 			     const int) = 0;

  virtual PdbCalBank* fetchBank(const std::string &,
				PdbBankID,
				const std::string &,
				const PHTimeStamp &) = 0;

  virtual PdbCalBank* fetchClosestBank(const std::string &,
				       PdbBankID,
				       const std::string &,
				       PHTimeStamp &) = 0;  

  // virtual void fetchAllBanks(PdbBankList &,
  // 			     const std::string &,
  // 			     PdbBankID,
  // 			     const std::string &,
  // 			     PHTimeStamp &) = 0;

  // virtual void fetchAllBanks(PdbBankList &,
  // 			     const std::string &,
  // 			     const std::string &,
  // 			     PHTimeStamp &) = 0;

  virtual PdbApplication* getApplication() = 0;

  virtual void fillCalibObject(PdbCalBank*,
			       const std::string &,
			       PHTimeStamp &) = 0;

  virtual void GetUsedBankRids(std::map<std::string,std::set<int> > &/*usedbanks*/) const {}
  virtual void ClearUsedBankRids() {}
  virtual  void SetMaxInsertTime(const PHTimeStamp &/*tMax*/) {}

protected:

  static  PdbBankManager *__instance; 

};

#endif /* PDBCAL_BASE_PDBBANKMANAGER_H */

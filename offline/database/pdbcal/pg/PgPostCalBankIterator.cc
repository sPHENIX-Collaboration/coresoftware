#include "PgPostCalBankIterator.h"
#include "PgPostApplication.h"
#include "PgPostBankManager.h"
#include "PgPostBankWrapper.h"

#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbCalBank.h>

#include <phool/PHTimeStamp.h>          // for PHTimeStamp, operator<<
#include <phool/phool.h>

#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLStatement.h>

#include <TString.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>


using namespace std;

//_____________________________________________________________________________
PgPostCalBankIterator::PgPostCalBankIterator(PgPostBankManager& bm)
  : fBM(bm)
  , fApplication(0)
  , fDBName("")
  , fTableName("")
  , fIsValid(false)
  , fBankID(-1)
  , fSQLStatement(0)
  , fSQLResultSet(0)
{
  fApplication = dynamic_cast<PgPostApplication*>(fBM.getApplication());
  if (!fApplication)
  {
    cout << PHWHERE << "dynamic_cast failed, exiting" << endl;
    exit(1);
  }
}

//_____________________________________________________________________________
PgPostCalBankIterator::~PgPostCalBankIterator()
{
  delete fSQLResultSet;
  delete fSQLStatement;
}

//_____________________________________________________________________________
bool PgPostCalBankIterator::init(const string& fulldbname, const PdbBankID& bankid)
{
  fDBName = fulldbname;
  fBankID = bankid;

  if (!fApplication->startRead())
  {
    std::cerr << PHWHERE << "Cannot start a read transaction"
              << std::endl;
    return false;
  }

  TSQLConnection* con = fApplication->getConnection();
  if (!con)
  {
    std::cout << PHWHERE << " Cannot get TSQLConnection, exiting" << std::endl;
    exit(1);
  }

  fSQLStatement = con->CreateStatement();

  std::ostringstream query;

  query << "select * from " << fulldbname;

  if (fBankID.getInternalValue() != -1)
  {
    query << " where bankID = "
          << bankid.getInternalValue();
  }

  fSQLResultSet = fSQLStatement->ExecuteQuery(query.str().c_str());

  fIsValid = false;

  if (fSQLResultSet)
  {
    fIsValid = true;
  }

  return fIsValid;
}

//_____________________________________________________________________________
bool PgPostCalBankIterator::isValid() const
{
  return fIsValid;
}

//_____________________________________________________________________________
PdbCalBank*
PgPostCalBankIterator::next()
{
  if (!fIsValid)
  {
    return 0;
  }

  //
  // Search the complete database for the next banks with
  // matching validity range(s).
  //

  ValPeriod insert(0, 0);
  ValPeriod end(0, 0);
  ValPeriod start(0, 0);

  const TimeMap::iterator it = fTimeMap.find("InsertTime");

  bool insertIntervalGiven = (it != fTimeMap.end());
  if (insertIntervalGiven)
  {
    insert = it->second;
  }

  const TimeMap::iterator it2 = fTimeMap.find("EndVal");
  bool endIntervalGiven = (it2 != fTimeMap.end());
  if (endIntervalGiven)
  {
    end = it2->second;
  }

  const TimeMap::iterator it3 = fTimeMap.find("StartVal");
  bool startIntervalGiven = (it3 != fTimeMap.end());
  if (startIntervalGiven)
  {
    start = it3->second;
  }

  while (fSQLResultSet->Next())
  {
    if (insertIntervalGiven)
    {
      // skip banks not corresponding to required insert interval
      time_t inserttime = fSQLResultSet->GetLong(2);
      if (inserttime >= insert.end() ||
          inserttime < insert.start())
      {
        continue;
      }
    }
    if (startIntervalGiven)
    {
      // skip banks not corresponding to required start interval
      time_t starttime = fSQLResultSet->GetLong(3);
      if (starttime > start.end() ||
          starttime < start.start())
      {
        continue;
      }
    }

    if (endIntervalGiven)
    {
      // skip banks not corresponding to required end interval
      time_t endtime = fSQLResultSet->GetLong(4);
      if (endtime > end.end() ||
          endtime < end.start())
      {
        continue;
      }
    }

    // Found a good match
    PdbCalBank* bank = static_cast<PdbCalBank*>(fSQLResultSet->GetObject(7));
    PgPostBankWrapper* bw = new PgPostBankWrapper(bank);
    bw->setBankID(fSQLResultSet->GetInt(1));
    bw->setInsertTime(fSQLResultSet->GetLong(2));
    bw->setStartValTime(fSQLResultSet->GetLong(3));
    bw->setEndValTime(fSQLResultSet->GetLong(4));
    bw->setDescription(fSQLResultSet->GetString(5).Data());
    bw->setUserName(fSQLResultSet->GetString(6).Data());
    bw->setTableName(fTableName);
    return bw;
  }

  fIsValid = false;
  return 0;
}

//_____________________________________________________________________________
void PgPostCalBankIterator::print(std::ostream& os) const
{
  TimeMap::const_iterator it;

  os << "PgPostCalBankIterator for database "
     << fDBName << " and BankID " << fBankID.getInternalValue()
     << " (table=" << fTableName << ")"
     << std::endl;

  for (it = fTimeMap.begin(); it != fTimeMap.end(); ++it)
  {
    os << it->first << "=[" << PHTimeStamp(it->second.start())
       << "," << PHTimeStamp(it->second.end()) << "]" << std::endl;
  }
}

//_____________________________________________________________________________
void PgPostCalBankIterator::setBankID(const PdbBankID& id)
{
  fBankID = id;
}

//_____________________________________________________________________________
void PgPostCalBankIterator::setEndValTimeLimits(const PHTimeStamp& min,
                                                const PHTimeStamp& max)
{
  fTimeMap["EndVal"] = ValPeriod(min, max);
}

//_____________________________________________________________________________
void PgPostCalBankIterator::setInsertTimeLimits(const PHTimeStamp& min,
                                                const PHTimeStamp& max)
{
  fTimeMap["InsertTime"] = ValPeriod(min, max);
}

//_____________________________________________________________________________
void PgPostCalBankIterator::setStartValTimeLimits(const PHTimeStamp& min,
                                                  const PHTimeStamp& max)
{
  fTimeMap["StartVal"] = ValPeriod(min, max);
}

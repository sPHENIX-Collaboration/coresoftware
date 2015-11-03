#include "PgPostBankManager.hh"
#include "PgPostHelper.hh"
#include "PgPostCalBankIterator.hh"
#include "PgPostBankWrapper.hh"
#include "PgPostApplication.hh"
#include "PgPostCalBank.hh"
#include "RunToTimePg.hh"

#include <pdbcalbase/PdbBankID.hh>
#include <pdbcalbase/PdbBankList.hh>
#include <pdbcalbase/PdbBankManagerFactory.hh>
#include <pdbcalbase/PdbCalBank.hh>
#include <pdbcalbase/PdbClassMap.hh>

#include <phool/PHPointerList.h>

#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLResultSetMetaData.h>
#include <RDBC/TSQLPreparedStatement.h>
#include <RDBC/TSQLDatabaseMetaData.h>

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>

using namespace std;

namespace
{

  PdbBankManager* singletonCreator()
  {

    // rememeber that this will not neccessarily return a
    // pointer to the singleton PgPostBankManager. If 
    // an Objy application is instantiated, it will return 0.
    return PgPostBankManager::instance();
  }

  const std::string name = "Pg";
  const bool registered =
    PdbBankManagerFactory::instance().registerCreator(name, singletonCreator, "PdbBankManager");
}

PgPostBankManager *PgPostBankManager::mySpecificCopy = 0;

PgPostBankManager* PgPostBankManager::instance()
{
  return mySpecificCopy;

}

int PgPostBankManager::Register()
{
  if ( __instance ) return -1;
  mySpecificCopy  = new  PgPostBankManager();
  __instance = mySpecificCopy;
  return 0;
}



PgPostBankManager::PgPostBankManager()
{
  tMaxInsertTime.setToFarFuture();
#ifdef DEBUG
  cout << PHWHERE << endl;
#endif
}

PgPostBankManager::~PgPostBankManager()
{
  mySpecificCopy = 0;
}

PdbCalBankIterator*
PgPostBankManager::getIterator()
{
  return new PgPostCalBankIterator(*this);
}

PdbCalBank* PgPostBankManager::createBank(const int beginRunNumber, const int endRunNumber, const char *className, PdbBankID bankID, const char *description, const char *bankName)
{
  RunToTime *runTime = RunToTime::instance();

  PHTimeStamp *beginRunTime = runTime->getBeginTime(beginRunNumber);
  if (beginRunTime != 0)
    {
      PHTimeStamp startTime = *(beginRunTime);

      PHTimeStamp *endRunTime = runTime->getEndTime(endRunNumber);
      if (endRunTime != 0)
	{
	  PHTimeStamp endTime = *(endRunTime);
	  delete beginRunTime;
	  delete endRunTime;
	  if (startTime >= endTime)
	    {
	      cout << PHWHERE << "Bad Start/EndRun Time: Start Time: "
		   << startTime << " >= End Time: "
		   << endTime << endl;
	      return 0;
	    }
	  return createBank(className, bankID, description, startTime, endTime, bankName);
	}
      else
	{
          delete beginRunTime;
	  cout << PHWHERE << "endTime = 0" << endl;
	  return 0;
	}
      delete beginRunTime;
    }
  else
    {
      cout << PHWHERE << "beginTime = 0" << endl;
    }
  return 0;
}


PdbCalBank* PgPostBankManager::createBank(const int runNumber, const char *className, PdbBankID bankID, const char *description, const char *bankName, const time_t duration)
{
  RunToTime *runTime = RunToTime::instance();

  PHTimeStamp *runBeginTime = runTime->getBeginTime(runNumber);
  if (runBeginTime != 0)
    {
      PHTimeStamp startTime = *(runBeginTime);
      PHTimeStamp endTime = startTime;

      if (duration == 0)
        {
          // duration == 0 is "flag" for end of unix time!!
          endTime = PHTimeStamp(2038, 1, 17, 0, 0, 0);
        }
      else
        {
          endTime += duration;
        }
      delete runBeginTime;
      return createBank(className, bankID, description, startTime, endTime, bankName);
    }
  return 0;
}


PdbCalBank* PgPostBankManager::createBank(const char *className, PdbBankID bankID, const char* descr, PHTimeStamp & tStart, PHTimeStamp & tStop, const char *tablename)
{

  string realName = getRealName(className);
  const char *rName = realName.c_str();
  PdbClassMap<PdbCalBank> *classMap = PdbClassMap<PdbCalBank>::instance();
  if (classMap->find(rName) != classMap->end())
    {
      std::string a = getTableName(tablename);
      PdbCalBank *b = (*classMap)[rName];
      PdbCalBank *b1 = b->clone();
      PgPostBankWrapper *bw = new PgPostBankWrapper(b1);
      bw->setBankID(bankID.getInternalValue());
      PHTimeStamp ts;
      bw->setInsertTime(ts);
      bw->setStartValTime(tStart);
      bw->setEndValTime(tStop);
      bw->setDescription(descr);
      bw->setUserName("pdbcal");
      bw->setTableName(a.c_str());
      return bw;
    }
  else
    {
      std::cerr << PHWHERE << "\t NO BANK " << rName 
		<< " IN THE MAP" << std::endl;
      return 0;
    }
}

PdbCalBank* PgPostBankManager::fetchBank(const char *className, PdbBankID bankID, const char *bankName, const int runNumber)
{
  RunToTime *runTime = RunToTime::instance();

  PHTimeStamp *runBeginTime = runTime->getBeginTime(runNumber);
  if (runBeginTime != 0)
    {
      PHTimeStamp searchTime = *(runBeginTime);
      delete runBeginTime;
      return fetchBank(className, bankID, bankName, searchTime);
    }
  return 0;
}


PdbCalBank* PgPostBankManager::fetchClosestBank(const char *className, PdbBankID bankID, const char *bankName, const int runNumber)
{
  RunToTime *runTime = RunToTime::instance();

  PHTimeStamp *runBeginTime = runTime->getBeginTime(runNumber);
  if (runBeginTime != 0)
    {
      PHTimeStamp searchTime = *(runBeginTime);
      delete runBeginTime;
      return fetchClosestBank(className, bankID, bankName, searchTime);
    }
  return 0;
}



// void PgPostBankManager::fetchAllBanks(PdbBankList & bankList, const char *className, PdbBankID bankID, const char *bankName, const int runNumber)
// {
//   RunToTime *runTime = RunToTime::instance();

//   PHTimeStamp *runBeginTime = runTime->getBeginTime(runNumber);
//   if (runBeginTime != 0)
//     {
//       PHTimeStamp searchTime = *(runBeginTime);
//       delete runBeginTime;
//       fetchAllBanks(bankList, className, bankID, bankName, searchTime);
//     }

// }


// void PgPostBankManager::fetchAllBanks(PdbBankList & bankList, const char *className, const char *bankName, const int runNumber)
// {
//   RunToTime *runTime = RunToTime::instance();

//   PHTimeStamp *runBeginTime = runTime->getBeginTime(runNumber);
//   if (runBeginTime != 0)
//     {
//       PHTimeStamp searchTime = *(runBeginTime);
//       delete runBeginTime;
//       fetchAllBanks(bankList, className, bankName, searchTime);
//     }
// }

//__________________________________________________________________________________
PdbCalBank* PgPostBankManager::fetchBank(const char *className, PdbBankID bankID, const char *bankName, const PHTimeStamp &searchTime)
{
#ifdef DEBUG
  cout << "Fetching " << className << " from " << bankName << endl;
#endif

  std::string a = getTableName(bankName);
  PgPostApplication* ap = PgPostApplication::instance();
  if (!ap)
    {
      cout << PHWHERE << " PgPostApplication instance is NULL, exiting" << endl;
      exit(1);
    }

  TSQLConnection *con = ap->getConnection();
  if (!con)
    {
      cout << PHWHERE << " Cannot get TSQLConnection, exiting" << endl;
      exit(1);
    }

  TSQLStatement *stmt = con->CreateStatement();
  time_t sT = searchTime.getTics();
  std::ostringstream tem;
  std::ostringstream t2;

  // cout << bankID.getInternalValue() << endl;
  t2 << "select * from " << a
     << " where bankID = " << bankID.getInternalValue()
     << " and startvaltime <= " << sT
     << " and endvaltime > " << sT;


  tem << "select * from ("
      << t2.str()
      << ") as foo where inserttime = "
      << "(select max(inserttime) from ("
      << t2.str()
      << " and inserttime <= "
      << tMaxInsertTime.getTics() 
      << ") as foobar)"
      << " order by rid desc";

#ifdef DEBUG
  cout << "exe : " << tem.str() << endl;
#endif

  std::auto_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));  
  if( (&*rs) && rs->Next())
  {
    PdbCalBank *bank = (PdbCalBank *)(rs->GetObject(7));
    PgPostBankWrapper* bw = new PgPostBankWrapper(bank);
    bw->setBankID(rs->GetInt(1));
    bw->setInsertTime(rs->GetLong(2));
    bw->setStartValTime(rs->GetLong(3));
    bw->setEndValTime(rs->GetLong(4));
    bw->setDescription(string(rs->GetString(5)));
    bw->setUserName(string(rs->GetString(6)));
    bw->setTableName(a);
    #ifdef DEBUG
    bw->printHeader();
    #endif
    int rid = rs->GetInt("rid");
    
    // insert new id in bank list matching name
    /*
    Remark: when the key "a" is not already in the map, it is inserted automatically by the call below,
    using the default constructor of the object associated to the key, here std::set<int> 
    */
    BankRid[a].insert(rid);
    return bw;
    
  } else {
    std::cerr << PHWHERE << "NO Bank found : " << tem.str() << std::endl;
    return 0;
  }
}






//__________________________________________________________________________________
PdbCalBank* PgPostBankManager::fetchClosestBank(const char *className, PdbBankID bankID, const char *bankName, PHTimeStamp &searchTime)
{ 
  cout << PHWHERE << " PdbBankManager::fetchClosestBank: This method is not implemented" << endl;
  exit(1);
}



// void PgPostBankManager::fetchAllBanks(PdbBankList & bankList, const char *className, PdbBankID bankID, const char *bankName, PHTimeStamp &searchTime)
// {
// #ifdef DEBUG
//   cout << "Fetching " << className << " from " << bankName << endl;
// #endif

//   std::string tablename = getTableName(bankName);

//   TSQLResultSet *rs;
//   PgPostBankWrapper *bw;

//   PgPostApplication *ap = PgPostApplication::instance();
//   if (!ap)
//     {
//       cout << PHWHERE << " PgPostApplication instance is NULL, exiting" << endl;
//       exit(1);
//     }
//   TSQLConnection *con = ap->getConnection();
//   if (!con)
//     {
//       cout << PHWHERE << " Cannot get TSQLConnection, exiting" << endl;
//       exit(1);
//     }
//   TSQLStatement* stmt = con->CreateStatement();
//   if (!stmt)
//     {
//       cout << PHWHERE << " Cannot create TSQL statement, exiting" << endl;
//       exit(1);
//     }

//   time_t sT = searchTime.getTics();
//   ostringstream tem;
//   tem << "select * from "  << tablename << " where bankID = "
//       << bankID.getInternalValue() << " and startvaltime <= "
//       << sT << " and endvaltime > " << sT;
//   cout << "exe: " << tem.str() << endl;
//   rs = stmt->ExecuteQuery(tem.str());
//   while (rs->Next())
//     {
//       bw = static_cast<PgPostBankWrapper*>(rs->GetObject(7));
// #ifdef DEBUG
//       bw->printHeader();
// #endif
//       bankList.append(bw);
//     }
//   cout << "bankList len= " << bankList.length() << endl;
// }




// void PgPostBankManager::fetchAllBanks(PdbBankList & bankList, const char *className, const char *bankName, PHTimeStamp &searchTime)
// {
// #ifdef DEBUG
//   cout << "Fetching " << className << " from " << bankName << endl;
// #endif

//   std::string tablename = getTableName(bankName);

//   TSQLResultSet *rs;
//   PgPostBankWrapper *bw;

//   PgPostApplication *ap = PgPostApplication::instance();
//   if (!ap)
//     {
//       cout << PHWHERE << " PgPostApplication instance is NULL, exiting" << endl;
//       exit(1);
//     }

//   TSQLConnection *con = ap->getConnection();
//   if (!con)
//     {
//       cout << PHWHERE << " Cannot get TSQLConnection, exiting" << endl;
//       exit(1);
//     }

//   TSQLStatement* stmt = con->CreateStatement();
//   if (!stmt)
//     {
//       cout << PHWHERE << " Cannot create TSQL statement, exiting" << endl;
//       exit(1);
//     }

//   time_t sT = searchTime.getTics();
//   ostringstream tem;
//   tem << "select * from " << tablename << " where startvaltime <= "
//       << sT << " and endvaltime > " << sT;
//   cout << "exe: " << tem.str() << endl;
//   rs = stmt->ExecuteQuery(tem.str());
//   while (rs->Next())
//     {
//       bw = static_cast<PgPostBankWrapper*>(rs->GetObject(7));
// #ifdef DEBUG
//       bw->printHeader();
// #endif
//       bankList.append(bw);
//     }
//   cout << "bankList len= " << bankList.length() << endl;
// }





PdbApplication* PgPostBankManager::getApplication(PHBoolean pJob)
{
  return PgPostApplication::instance();
}

string 
PgPostBankManager::getRealName(const string &searchName)
{
  string realName = searchName;
  string pdbsubstring = "Pdb";
  realName.replace(realName.find(pdbsubstring),realName.find(pdbsubstring)+pdbsubstring.size(),"PgPost");
  return realName;
}

void
PgPostBankManager::GetUsedBankRids(map<string,set<int> > &usedbanks) const
{
  usedbanks = BankRid;
//   map<string,set<int> >::const_iterator bankiter;
//   for (bankiter = usedbanks.begin(); bankiter != usedbanks.end(); bankiter++)
//     {
//       cout << "BankName: " << bankiter->first << endl;
//       set<int>:: const_iterator siter;
//       for (siter = (bankiter->second).begin(); siter != (bankiter->second).end(); siter++)
// 	{
// 	  cout << "rid: " << *siter << endl;
// 	}
//     }
  return;
}

void 
PgPostBankManager::SetMaxInsertTime(const PHTimeStamp &tMax)
{
  cout << "Setting latest inserttime for calibrations to " << tMax 
       << " (" << tMax.getTics() << ")" << endl;
  tMaxInsertTime = tMax;
 return;
}

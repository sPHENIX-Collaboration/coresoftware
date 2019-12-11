#include "RunToTimePg.h"

#include <pdbcalbase/RunToTime.h>  // for RunToTime::__instance

#include <phool/PHTimeStamp.h>
#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>
#include <odbc++/types.h>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>                // for reverse_iterator
#include <sstream>
#include <string>
#include <utility>                 // for pair

using namespace odbc;
using namespace std;

// cache entries for 10 runs which should be sufficient for everybody
unsigned int maxentries = 10;

RunToTimePg* RunToTimePg::mySpecificCopy = nullptr;

RunToTimePg::RunToTimePg()
{
  con = nullptr;
  return;
}

RunToTimePg::~RunToTimePg()
{
  mySpecificCopy = nullptr;
  __instance = nullptr;
  while (beginruntimes.begin() != beginruntimes.end())
  {
    delete beginruntimes.begin()->second;
    beginruntimes.erase(beginruntimes.begin());
  }
  while (endruntimes.begin() != endruntimes.end())
  {
    delete endruntimes.begin()->second;
    endruntimes.erase(endruntimes.begin());
  }
  return;
}

int RunToTimePg::Register()
{
  if (__instance)
  {
    return -1;
  }
  mySpecificCopy = new RunToTimePg();
  __instance = mySpecificCopy;
  return 0;
}

int RunToTimePg::GetConnection()
{
  if (!con)
  {
    try
    {
      con = DriverManager::getConnection("daq", "phnxrc", "");
    }
    catch (SQLException& e)
    {
      cout << PHWHERE
           << " Fatal Exception caught during DriverManager::getConnection" << endl;
      cout << "Message: " << e.getMessage() << endl;
      exit(1);
    }
  }
  return 0;
}

int RunToTimePg::DisconnectDB()
{
  delete con;
  con = nullptr;
  return 0;
}

PHTimeStamp*
RunToTimePg::getTime(const int runNumber, const string& what)
{
  PHTimeStamp* whatTime = nullptr;
  //  Establish connection to Postgres...
  GetConnection();  // on error this method will exit
  Statement* stmt = con->createStatement();

  std::ostringstream cmd;
  cmd << "select brunixtime,erunixtime,updateunixtime from run where runnumber = " << runNumber
      << " limit 1";

#ifdef DEBUG

  cout << cmd.str() << endl;
#endif

  //  Get results of search...
  ResultSet* rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException& e)
  {
    cout << "Fatal Exception caught during stmt->executeQuery(" << cmd.str() << ")" << endl;
    cout << "Message: " << e.getMessage() << endl;
    exit(1);
  }

  //  Fill stl maps with timestamps
  if (rs->next())
  {
    // if the maps grow too large we need to ditch them (e.g. a process in the online monitoring
    // which runs over a lot of runs. Chances are the last entry is the one for the current run
    // and this one might be reused. Chances are the lower run(s) have some special meaning
    // (forcing the use of a special calibration) and might be reused as well. So we remove the second
    // largest runnumber in the list
    if (beginruntimes.size() > maxentries)
    {
      // this picks the last run
      map<const int, PHTimeStamp*>::reverse_iterator iter = beginruntimes.rbegin();
      // go to the run before
      ++iter;
      int delrun = iter->first;
      delete iter->second;
      beginruntimes.erase(delrun);
      map<const int, PHTimeStamp*>::iterator iter2;
      iter2 = endruntimes.find(delrun);
      if (iter2 != endruntimes.end())
      {
        delete iter2->second;
        endruntimes.erase(iter2);
      }
    }
    // begin run time
    try
    {
      whatTime = new PHTimeStamp(rs->getInt("brunixtime"));
    }
    catch (SQLException& e)
    {
      cout << "Fatal Exception caught during rs->getInt(\"brunixtime\")" << endl;
      cout << "Message: " << e.getMessage() << endl;
      exit(1);
    }
    beginruntimes[runNumber] = whatTime;
    // end run time
    try
    {
      unsigned int eruntics = rs->getInt("erunixtime");
      if (eruntics == 0)

      {
        eruntics = rs->getInt("updateunixtime");
      }
      whatTime = new PHTimeStamp(eruntics);
    }
    catch (SQLException& e)
    {
      cout << "Fatal Exception caught during rs->getInt(\"brunixtime\")" << endl;
      cout << "Message: " << e.getMessage() << endl;
      exit(1);
    }
    endruntimes[runNumber] = whatTime;
    if (what == "brunixtime")
    {
      delete rs;
      return beginruntimes[runNumber];
    }
    else if (what == "erunixtime")
    {
      delete rs;
      return endruntimes[runNumber];
    }
    else
    {
      cout << "invalid time selection " << what << endl;
      exit(1);
    }
  }
  delete rs;
  return 0;
}

PHTimeStamp*
RunToTimePg::getBeginTime(const int runNumber)
{
  PHTimeStamp* BeginRunTime;
  map<const int, PHTimeStamp*>::const_iterator iter = beginruntimes.find(runNumber);
  if (iter == beginruntimes.end())
  {
    BeginRunTime = getTime(runNumber, "brunixtime");
  }
  else
  {
    BeginRunTime = iter->second;
  }
  if (!BeginRunTime)
  {
    return nullptr;
  }
  PHTimeStamp* TS = new PHTimeStamp(*BeginRunTime);
  return TS;
}

PHTimeStamp*
RunToTimePg::getEndTime(const int runNumber)
{
  PHTimeStamp* EndRunTime;
  map<const int, PHTimeStamp*>::const_iterator iter = endruntimes.find(runNumber);
  if (iter == endruntimes.end())
  {
    EndRunTime = getTime(runNumber, "erunixtime");
  }
  else
  {
    EndRunTime = iter->second;
  }
  if (!EndRunTime)
  {
    return nullptr;
  }
  PHTimeStamp* TS = new PHTimeStamp(*EndRunTime);
  return TS;
}

int RunToTimePg::getRunNumber(const PHTimeStamp& ts)
{
  GetConnection();
  Statement* stmt = con->createStatement();

  //  Make search string...
  std::ostringstream cmd;

  Timestamp timestp(ts.getTics());

  time_t tics = ts.getTics();

  cmd << "select runnumber from run where brunixtime <= " << tics
      << " order by runnumber desc limit 1";

  //  Get results of search...
  ResultSet* rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str().c_str());
  }
  catch (SQLException& e)
  {
    cout << "Fatal Exception caught during stmt->executeQuery(" << cmd.str() << ")" << endl;
    cout << "Message: " << e.getMessage() << endl;
    exit(1);
  }

  int runnumber = -1;
  time_t eruntics = 0;
  //  Set runnumber based upon return from
  if (rs->next())
  {
    runnumber = rs->getInt("runnumber");
  }
  delete rs;
  if (runnumber == -1)
  {
    goto cleanup;
  }

  cmd.str("");
  cmd << "select erunixtime from run where runnumber = " << runnumber;

  try
  {
    rs = stmt->executeQuery(cmd.str().c_str());
  }
  catch (SQLException& e)
  {
    cout << "Fatal Exception caught during stmt->executeQuery(" << cmd.str() << ")" << endl;
    cout << "Message: " << e.getMessage() << endl;
    exit(1);
  }

  //  Set runnumber based upon return from
  if (rs->next())
  {
    eruntics = rs->getInt("erunixtime");
    if (eruntics > 0 && eruntics < tics)
    {
      cout << "Timestamp " << tics
           << " not covered by any run, closest smaller begin run time is from run "
           << runnumber << endl;
      runnumber = -1;
    }
  }
  delete rs;

cleanup:
  return runnumber;
}

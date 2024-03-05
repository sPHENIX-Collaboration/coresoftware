#include "OnCalServer.h"
#include "OnCal.h"
#include "OnCalDBCodes.h"
#include "OnCalHistoBinDefs.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <phool/phool.h>

#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankList.h>
#include <pdbcalbase/PdbBankListIterator.h>
#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pgcal/PgPostCalBank.h>

#include <pdbcalbase/RunToTime.h>

#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>

// odbc++ classes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/preparedstatement.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#pragma GCC diagnostic pop

#include <boost/foreach.hpp>

#include <sys/utsname.h>
#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <sstream>

using namespace std;
using namespace odbc;

const static string cvstag = "OnCalv86";

odbc::Connection *DBconnection = nullptr;

OnCalServer *OnCalServer::instance()
{
  if (__instance)
  {
    OnCalServer *oncal = dynamic_cast<OnCalServer *>(__instance);
    return oncal;
  }
  __instance = new OnCalServer();
  OnCalServer *oncal = dynamic_cast<OnCalServer *>(__instance);
  return oncal;
}

//---------------------------------------------------------------------

OnCalServer::OnCalServer(const std::string &name)
  : Fun4AllServer(name)
  , testmode(0)
  , recordDB(false)
  , SetEndTimeStampByHand(false)
  , SetBeginTimeStampByHand(false)
  , successTable("OnCal_status")
  , runNum(0)
  , nEvents(0)
  , eventcheckfrequency(1000)
  , database("calBookKeep")
{
  beginTimeStamp.setTics(0);
  endTimeStamp.setTics(0);

  OnCalServerVars = new TH1D("OnCalServerVars", "OnCalServerVars", OnCalHistoBinDefs::LASTBINPLUSONE, -0.5, (OnCalHistoBinDefs::LASTBINPLUSONE - 0.5));
  Fun4AllServer::registerHisto(OnCalServerVars);
  return;
}
//---------------------------------------------------------------------

OnCalServer::~OnCalServer()
{
  delete DBconnection;
  return;
}
//---------------------------------------------------------------------

PHTimeStamp *
OnCalServer::GetEndValidityTS()
{
  if (endTimeStamp.getTics())
  {
    PHTimeStamp *ts = new PHTimeStamp(endTimeStamp);
    return ts;
  }
  else
  {
    cout << PHWHERE << "Screwup - the end validity time is not set" << endl;
    exit(1);
  }
}
//---------------------------------------------------------------------

PHTimeStamp *OnCalServer::GetBeginValidityTS()
{
  if (beginTimeStamp.getTics())
  {
    PHTimeStamp *ts = new PHTimeStamp(beginTimeStamp);
    return ts;
  }
  else
  {
    cout << PHWHERE << "Screwup - the begin validity time is not set" << endl;
    exit(1);
  }
}
//---------------------------------------------------------------------

void OnCalServer::dumpHistos()
{
  ostringstream filename;
  string fileprefix = "./";

  if (getenv("ONCAL_SAVEDIR"))
  {
    fileprefix = getenv("ONCAL_SAVEDIR");
    fileprefix += "/";
  }

  int compress = 3;
  map<string, set<string> >::const_iterator iter;
  //  map<string, TH1 *>::const_iterator hiter;
  TH1 *histo;
  set<string>::const_iterator siter;
  for (iter = calibratorhistomap.begin(); iter != calibratorhistomap.end(); ++iter)
  {
    filename.str("");
    filename << fileprefix << "Run_"
             << RunNumber()
             << "_" << iter->first << ".root";
    TFile *hfile = new TFile(filename.str().c_str(), "RECREATE",
                             "Created by Online Calibrator", compress);
    cout << "OnCalServer::dumpHistos() Output root file: " << filename.str() << endl;
    for (siter = (iter->second).begin(); siter != (iter->second).end(); ++siter)
    {
      histo = dynamic_cast<TH1 *>(getHisto(*siter));
      if (histo)
      {
        histo->Write();
      }
      else
      {
        cout << PHWHERE << "Histogram "
             << *siter << " not found, will not be saved in "
             << filename.str() << endl;
      }
    }
    hfile->Close();

    delete hfile;
  }
  return;
}

void OnCalServer::registerHisto(TH1 *h1d, OnCal *Calibrator, const int replace)
{
  if (Calibrator)
  {
    string calibratorname = Calibrator->Name();
    map<string, set<string> >::iterator iter;
    iter = calibratorhistomap.find(calibratorname);
    if (iter != calibratorhistomap.end())
    {
      (iter->second).insert(h1d->GetName());
    }
    else
    {
      set<string> newset;
      newset.insert(h1d->GetName());
      newset.insert("OnCalServerVars");
      calibratorhistomap[calibratorname] = newset;
    }
  }
  Fun4AllServer::registerHisto(h1d, replace);
  return;
}

void OnCalServer::unregisterHisto(const string &calibratorname)
{
  calibratorhistomap.erase(calibratorname);
  return;
}

int OnCalServer::process_event()
{
  Fun4AllServer::process_event();
  int i = 0;
  nEvents++;
  if ((nEvents % eventcheckfrequency) == 0)  // check every 1000 events
  {
    cout << nEvents << " events, testing" << endl;
    unsigned int j = 0;
    unsigned int ical = 0;
    vector<pair<SubsysReco *, PHCompositeNode *> >::const_iterator iter;
    for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
    {
      OnCal *oncal = dynamic_cast<OnCal *>(iter->first);
      if (oncal)
      {
        ical++;
        cout << "Name: " << oncal->Name()
             << " is " << oncal->AllDone() << endl;
        j += oncal->AllDone();
      }
    }
    if (j == ical)
    {
      cout << "Everyone is done after "
           << nEvents << " Events" << endl;
      i = 1;
    }
  }
  return i;
}

int OnCalServer::BeginRun(const int runno)
{
  if (runno <= 0)
  {
    cout << PHWHERE << "Invalid Run Number: " << runno << endl;
    exit(1);
  }
  FillRunListFromFileList();
  recoConsts *rc = recoConsts::instance();
  // we stick to the first runnumber, but after inheriting from
  // Fun4All we get a EndRun/BeginRun when the run number changes
  // so we have to catch this here
  if (RunNumber() != 0)
  {
    rc->set_IntFlag("RUNNUMBER", RunNumber());  // set rc flag back to previous run
    analysed_runs.push_back(runno);
    return 0;
  }
  RunNumber(runno);
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  // copy the subsys reco pointers to another set for
  // easier search (we only need the pointers to find
  // the subsystems with special timestamp/runnumber needs
  set<SubsysReco *> NeedOtherTimeStamp;
  map<string, set<SubsysReco *> >::const_iterator miter;
  set<SubsysReco *>::const_iterator siter;
  for (miter = requiredCalibrators.begin();
       miter != requiredCalibrators.end(); ++miter)
  {
    for (siter = miter->second.begin(); siter != miter->second.end(); ++siter)
    {
      NeedOtherTimeStamp.insert(*siter);
    }
  }

  int iret;
  int i = 0;
  int oncalrun = runno;
  int fun4allrun = runno;

  RunToTime *runTime = RunToTime::instance();
  PHTimeStamp *ts = runTime->getBeginTime(fun4allrun);
  PHTimeStamp OnCalBORTimeStamp = *ts;
  PHTimeStamp Fun4AllBORTimeStamp(OnCalBORTimeStamp);
  delete ts;
  if (requiredCalibrators.size() > 0)
  {
    fun4allrun = FindClosestCalibratedRun(runno);
    ts = runTime->getBeginTime(fun4allrun);
    Fun4AllBORTimeStamp = *ts;
    delete ts;
  }

  // we have to do the same TDirectory games as in the Init methods
  // save the current dir, cd to the subsystem name dir (which was
  // created in init) call the InitRun of the module and cd back

  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();
  set<string> droplist;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    OnCal *oncal = dynamic_cast<OnCal *>((*iter).first);
    if (oncal)
    {
      string table = "OnCal";
      table += (*iter).first->Name();
      check_create_subsystable(table);
      insertRunNumInDB(table, runNum);
      string calibname = (*iter).first->Name();
      add_calibrator_to_statustable(calibname);
      set<int>::const_iterator runiter;
      int calibstatus = GetCalibStatus(calibname, runNum);
      if (calibstatus > 0 && testmode == 0)
      {
        cout << calibname << " already ran for run " << runNum << endl;
        droplist.insert(calibname);
        unregisterSubsystem(oncal);
        unregisterHisto(calibname);
      }
      else
      {
        ostringstream stringarg;
        stringarg << OnCalDBCodes::STARTED;
        for (runiter = runlist.begin(); runiter != runlist.end(); ++runiter)
        {
          updateDB(successTable, calibname, stringarg.str(), *runiter);
        }
      }
    }
    if (NeedOtherTimeStamp.find((*iter).first) != NeedOtherTimeStamp.end())
    {
      cout << "changing timestamp for " << (*iter).first->Name() << endl;
      rc->set_IntFlag("RUNNUMBER", fun4allrun);
      // rc->set_TimeStamp(Fun4AllBORTimeStamp);
    }
    else
    {
      rc->set_IntFlag("RUNNUMBER", oncalrun);
      // rc->set_TimeStamp(OnCalBORTimeStamp);
    }
    if (droplist.find((*iter).first->Name()) == droplist.end())
    {
      iret = (*iter).first->InitRun(TopNode);
      if (iret == Fun4AllReturnCodes::ABORTRUN)
      {
        cout << PHWHERE << "Module " << (*iter).first->Name() << " issued Abort Run, exiting" << endl;
        exit(-1);
      }
      i += iret;
    }
  }
  gROOT->cd(currdir.c_str());

  rc->set_IntFlag("RUNNUMBER", oncalrun);
  //  rc->set_TimeStamp(OnCalBORTimeStamp);
  if (OnCalServerVars->GetBinContent(OnCalHistoBinDefs::FIRSTRUNBIN) == 0)
  {
    OnCalServerVars->SetBinContent(OnCalHistoBinDefs::FIRSTRUNBIN, runno);
    OnCalServerVars->SetBinContent(OnCalHistoBinDefs::BORTIMEBIN, (Stat_t) OnCalBORTimeStamp.getTics());
  }
  OnCalServerVars->SetBinContent(OnCalHistoBinDefs::LASTRUNBIN, (Stat_t) runno);
  ts = runTime->getEndTime(runno);
  if (ts)
  {
    OnCalServerVars->SetBinContent(OnCalHistoBinDefs::EORTIMEBIN, (Stat_t) ts->getTics());
    delete ts;
  }

  // disconnect from DB to save resources on DB machine
  // PdbCal leaves the DB connection open (PdbCal will reconnect without
  // problem if neccessary)
  DisconnectDB();
  // finally drop calibrators which have run already from module list
  unregisterSubsystemsNow();
  return i;
}

int OnCalServer::End()
{
  if (nEvents == 0)
  {
    cout << "No Events read, you probably gave me an empty filelist" << endl;
    return -1;
  }
  PdbApplication *app = PdbApplication::instance();
  app->setDBName("oncal");
  int i = 0;
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();

  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    else
    {
      if (Verbosity() > 2)
      {
        cout << "End: cded to " << newdirname.str().c_str() << endl;
      }
    }
    i += (*iter).first->End((*iter).second);
  }

  gROOT->cd(default_Tdirectory.c_str());
  currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    OnCal *oncal = dynamic_cast<OnCal *>((*iter).first);
    if (!oncal)
    {
      continue;
    }
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }

    string CalibratorName = oncal->Name();

    int verificationstatus = oncal->VerificationOK();
    int databasecommitstatus = oncal->CommitedToPdbCalOK();

    // report success database the status of the calibration
    if (recordDB)
    {
      string table = "OnCal";
      table += CalibratorName;

      ostringstream stringarg;
      if (databasecommitstatus == OnCalDBCodes::SUCCESS)
      {
        stringarg << OnCalDBCodes::COVERED;
      }
      else
      {
        stringarg << OnCalDBCodes::FAILED;
      }
      set<int>::const_iterator runiter;
      for (runiter = runlist.begin(); runiter != runlist.end(); ++runiter)
      {
        updateDB(successTable, CalibratorName, stringarg.str(), *runiter);
      }
      // update the first run which was used in the calibration
      // with the real status
      updateDB(successTable.c_str(), CalibratorName.c_str(),
               databasecommitstatus);

      stringarg.str("");
      stringarg << databasecommitstatus;
      updateDB(table, "committed", stringarg.str(), RunNumber());

      stringarg.str("");
      stringarg << verificationstatus;
      updateDB(table, "verified", stringarg.str(), RunNumber());

      odbc::Timestamp stp(time(nullptr));
      updateDB(table, "date", stp.toString(), RunNumber());
      updateDB(table, "comment", oncal->Comment(), RunNumber());
      time_t beginticks = beginTimeStamp.getTics();
      stringarg.str("");
      stringarg << beginticks;
      updateDB(table, "startvaltime", stringarg.str(), RunNumber());
      stp.setTime(beginticks);
      updateDB(table, "begintime", stp.toString(), RunNumber());
      time_t endticks = endTimeStamp.getTics();
      stringarg.str("");
      stringarg << endticks;
      updateDB(table, "endvaltime", stringarg.str(), RunNumber());
      stp.setTime(endticks);
      updateDB(table, "endtime", stp.toString(), RunNumber());

      string filelist = "";
      BOOST_FOREACH (Fun4AllSyncManager *sync, SyncManagers)
      {
        BOOST_FOREACH (Fun4AllInputManager *inmgr, sync->GetInputManagers())
        {
          BOOST_FOREACH (string infile, inmgr->GetFileOpenedList())
          {
            filelist += (infile).substr(((infile).find_last_of('/') + 1), (infile).size());
            filelist += " ";  // this needs to be stripped again for last entry
          }
        }
      }
      filelist.pop_back();  // strip empty space at end from loop
      cout << "FileList: " << filelist << endl;
      updateDB(table, "files", filelist, RunNumber());
      updateDB(table, "cvstag", cvstag, RunNumber());
    }

    cout << "SERVER SUMMARY: " << oncal->Name() << "  "
         << (verificationstatus == 1 ? "Verification: SUCCESS  " : "Verification: FAILURE  ")
         << (databasecommitstatus == 1 ? "DB commit: SUCCESS  " : "DB commit: FAILURE  ")
         << endl;

    printStamps();
  }
  gROOT->cd(currdir.c_str());
  dumpHistos();  // save the histograms in files
  return i;
}
//---------------------------------------------------------------------

void OnCalServer::Print(const string &what) const
{
  Fun4AllServer::Print(what);
  if (what == "ALL" || what == "CALIBRATOR")
  {
    // loop over the map and print out the content
    // (name and location in memory)

    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of Calibrators in OnCalServer:" << endl;

    vector<pair<SubsysReco *, PHCompositeNode *> >::const_iterator miter;
    for (miter = Subsystems.begin();
         miter != Subsystems.end(); ++miter)
    {
      OnCal *oncal = dynamic_cast<OnCal *>((*miter).first);
      if (oncal)
      {
        cout << oncal->Name() << endl;
      }
    }
    cout << endl;
  }
  if (what == "ALL" || what == "REQUIRED")
  {
    // loop over the map and print out the content
    // (name and location in memory)

    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of required Calibrations in OnCalServer:" << endl;

    map<string, set<SubsysReco *> >::const_iterator iter;
    set<SubsysReco *>::const_iterator siter;
    for (iter = requiredCalibrators.begin();
         iter != requiredCalibrators.end(); ++iter)
    {
      cout << iter->first << " calibrations are needed by " << endl;
      for (siter = iter->second.begin(); siter != iter->second.end(); ++siter)
      {
        cout << (*siter)->Name() << endl;
      }
    }
    cout << endl;
  }
  if (what == "ALL" || what == "FILES")
  {
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of PRDF Files in OnCalServer:" << endl;
    BOOST_FOREACH (Fun4AllSyncManager *sync, SyncManagers)
    {
      BOOST_FOREACH (Fun4AllInputManager *inmgr, sync->GetInputManagers())
      {
        BOOST_FOREACH (string infile, inmgr->GetFileList())
        {
          cout << "File: " << infile << endl;
        }
      }
    }
  }
  if (what == "ALL" || what == "RUNS")
  {
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of Run Numbers in OnCalServer:" << endl;
    set<int>::const_iterator liter;
    for (liter = runlist.begin(); liter != runlist.end(); ++liter)
    {
      cout << "Run : " << *liter << endl;
    }
  }

  return;
}

void OnCalServer::printStamps()
{
  cout << endl
       << endl;
  cout << "*******************************************" << endl;
  cout << "*    VALIDITY RANGE FOR THIS CALIBRATION  *" << endl;
  cout << "*                                         *" << endl;
  cout << "* Used Run     :   ";
  cout << runNum << endl;
  cout << endl;
  cout << "* Begin Valid  :   ";
  beginTimeStamp.print();
  cout << endl;
  cout << "* End Valid    :   ";
  endTimeStamp.print();
  cout << endl;
  cout << "*                                         *" << endl;
  cout << "*******************************************" << endl;
  cout << endl
       << endl
       << endl;
}

//---------------------------------------------------------------------

void OnCalServer::RunNumber(const int runnum)
{
  runNum = runnum;
  SetBorTime(runnum);
  if (recordDB)
  {
    set<int>::const_iterator runiter;
    time_t beginrunticks;
    time_t endrunticks;
    ostringstream stringarg;
    odbc::Timestamp stp;
    for (runiter = runlist.begin(); runiter != runlist.end(); ++runiter)
    {
      insertRunNumInDB(successTable, *runiter);
      GetRunTimeTicks(*runiter, beginrunticks, endrunticks);
      stringarg.str("");
      stringarg << beginrunticks;
      updateDB(successTable, "startvaltime", stringarg.str(), *runiter);
      stp.setTime(beginrunticks);
      updateDB(successTable, "beginrun", stp.toString(), *runiter);
      stringarg.str("");
      stringarg << endrunticks;
      updateDB(successTable, "endvaltime", stringarg.str(), *runiter);
      stp.setTime(endrunticks);
      updateDB(successTable, "endrun", stp.toString(), *runiter);
    }
  }
  if (runlist.size() > 0)
  {
    SetEorTime(*runlist.rbegin());
  }
  return;
}

//---------------------------------------------------------------------

bool OnCalServer::connectDB()
{
  if (DBconnection)
  {
    return true;
  }

  bool failure = true;
  int countdown = 10;
  while (failure && countdown > 0)
  {
    failure = false;
    try
    {
      DBconnection =
          DriverManager::getConnection(database.c_str(), "phnxrc", "");
    }
    catch (SQLException &e)
    {
      cout << "Cannot connect to " << database.c_str() << endl;
      cout << e.getMessage() << endl;
      cout << "countdown: " << countdown << endl;
      countdown--;
      failure = true;
      sleep(100);  // try again in 100 secs
    }
  }
  if (failure)
  {
    cout << "could not connect to DB after 10 tries in 1000 secs, giving up" << endl;
    exit(-1);
  }
  cout << "connected to " << database.c_str() << " database." << endl;
  return true;
}
//---------------------------------------------------------------------

int OnCalServer::DisconnectDB()
{
  delete DBconnection;
  DBconnection = nullptr;
  return 0;
}
//---------------------------------------------------------------------

bool OnCalServer::insertRunNumInDB(const string &DBtable, const int runno)
{
  if (findRunNumInDB(DBtable, runno))
  {
    return true;
  }

  cout << "new row will be created in DB for run " << runno << endl;

  Statement *statement = nullptr;
  statement = DBconnection->createStatement();
  ostringstream cmd;
  cmd << "INSERT INTO "
      << DBtable
      << " (runnumber) VALUES ("
      << runno << ")";

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::insertRunNumInDB() ... ";
    cout << "executing SQL statements ..." << endl;
    cout << cmd.str() << endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

bool OnCalServer::findRunNumInDB(const string &DBtable, const int runno)
{
  if (!DBconnection)
  {
    connectDB();
  }
  Statement *statement = nullptr;
  ResultSet *rs = nullptr;
  ostringstream cmd;
  cmd << "SELECT runnumber FROM "
      << DBtable
      << " WHERE runnumber = "
      << runno;

  statement = DBconnection->createStatement();

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::findRunNumInDB() ";
    cout << "executing SQL statement ..." << endl
         << cmd.str() << endl;
  }

  try
  {
    rs = statement->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << PHWHERE << " exception caught: " << e.getMessage() << endl;
    return false;
  }

  int entry = 0;
  if (rs->next())
  {
    try
    {
      entry = rs->getInt("runnumber");
    }
    catch (SQLException &e)
    {
      cout << PHWHERE << " exception caught: " << e.getMessage() << endl;
      return false;
    }
  }
  else
  {
    return false;
  }
  cout << "run number " << entry << " already exists in DB" << endl;
  return true;
}

bool OnCalServer::updateDBRunRange(const char *table, const char *column, const int entry, const int firstrun, const int lastrun)
{
  if (!DBconnection)
  {
    connectDB();
  }

  Statement *statement = nullptr;

  TString command = "UPDATE ";
  command += table;
  command += " SET ";
  command += column;
  command += " = ";
  command += entry;
  command += " WHERE runnumber >= ";
  command += firstrun;
  command += " and runnumber <= ";
  command += lastrun;

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::updateDB() ... ";
    cout << "executin SQL statement ... " << endl;
    cout << command.Data() << endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command.Data());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

bool OnCalServer::updateDB(const char *table, const char *column, int entry)
{
  if (!DBconnection)
  {
    connectDB();
  }

  Statement *statement = nullptr;

  TString command = "UPDATE ";
  command += table;
  command += " SET ";
  command += column;
  command += " = ";
  command += entry;
  command += " WHERE runnumber = ";
  command += runNum;

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::updateDB() ... ";
    cout << "executin SQL statement ... " << endl;
    cout << command.Data() << endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command.Data());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return false;
  }

  return true;
}
//---------------------------------------------------------------------

bool OnCalServer::updateDB(const char *table, const char *column, bool entry)
{
  if (!DBconnection)
  {
    connectDB();
  }
  Statement *statement = nullptr;

  TString command = "UPDATE ";
  command += table;
  command += " set ";
  command += column;
  command += " = '";
  command += static_cast<int>(entry);
  command += "' WHERE runnumber = ";
  command += runNum;

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::updateDB() ... ";
    cout << "executin SQL statement ... " << endl;
    cout << command.Data() << endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command.Data());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

int OnCalServer::updateDB(const string &table, const string &column,
                          const time_t ticks)
{
  if (!DBconnection)
  {
    connectDB();
  }
  Statement *statement = nullptr;

  ostringstream cmd;
  statement = DBconnection->createStatement();
  cmd << "UPDATE "
      << table
      << " set "
      << column
      << " = "
      << ticks
      << " WHERE runnumber = "
      << runNum;

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::updateDB() ... ";
    cout << "executin SQL statement ... " << endl;
    cout << cmd.str() << endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return -1;
  }
  return 0;
}
//---------------------------------------------------------------------

bool OnCalServer::updateDB(const string &table, const string &column,
                           const string &entry, const int runno, const bool append)
{
  if (!DBconnection)
  {
    connectDB();
  }

  Statement *statement = nullptr;

  statement = DBconnection->createStatement();

  string comment = "";
  ostringstream cmd;
  if (append)
  {
    ResultSet *rs = nullptr;
    ostringstream query;
    query << "SELECT * FROM "
          << table
          << " WHERE runnumber = "
          << runno;

    try
    {
      rs = statement->executeQuery(query.str());
    }
    catch (SQLException &e)
    {
      cout << "in function OnCalServer::updateDB() ... ";
      cout << "run number " << runno << "not found in DB" << endl;
      cout << e.getMessage() << endl;
    }

    rs->next();
    try
    {
      comment = rs->getString(column);
      comment += " ";  // add empty space between comments
    }
    catch (SQLException &e)
    {
      cout << "in function OnCalServer::updateDB() ... " << endl;
      cout << "nothing to append." << endl;
      cout << e.getMessage() << endl;
    }
    delete rs;
  }

  comment += entry;
  cmd << "UPDATE "
      << table
      << " set "
      << column
      << " = '"
      << comment
      << "' WHERE runnumber = "
      << runno;

  if (Verbosity() == 1)
  {
    cout << "in function OnCalServer::updateDB() ... ";
    cout << "executin SQL statement ... " << endl;
    cout << cmd.str() << endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << e.getMessage() << endl;
    return false;
  }
  delete statement;
  return true;
}

//---------------------------------------------------------------------

int OnCalServer::check_create_subsystable(const std::string &tablename)
{
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  vector<pair<string, string> > calibrator_columns;
  vector<pair<string, string> >::const_iterator coliter;
  calibrator_columns.emplace_back("runnumber", "int NOT NULL");
  calibrator_columns.emplace_back("verified", "int default -2");
  calibrator_columns.emplace_back("committed", "int default -2");
  calibrator_columns.emplace_back("date", "timestamp(0) with time zone");
  calibrator_columns.emplace_back("comment", "text");
  calibrator_columns.emplace_back("files", "text");
  calibrator_columns.emplace_back("cvstag", "text");
  calibrator_columns.emplace_back("startvaltime", "bigint");
  calibrator_columns.emplace_back("endvaltime", "bigint");
  calibrator_columns.emplace_back("begintime", "timestamp(0) with time zone");
  calibrator_columns.emplace_back("endtime", "timestamp(0) with time zone");

  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  cmd << "SELECT * FROM " << tablename << " LIMIT 1" << ends;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << tablename << " does not exist, will create it" << endl;
    //      cout << "Message: " << e.getMessage() << endl;
  }
  if (!rs)
  {
    cmd.str("");
    cmd << "CREATE TABLE "
        << tablename
        << "(";
    for (coliter = calibrator_columns.begin(); coliter != calibrator_columns.end(); ++coliter)
    {
      cmd << (*coliter).first << " " << (*coliter).second << ", ";
    }

    cmd << "primary key(runnumber))";
    stmt->executeUpdate(cmd.str());
  }
  else  // check if the all columns exist
  {
    for (coliter = calibrator_columns.begin(); coliter != calibrator_columns.end(); ++coliter)
    {
      try
      {
        rs->findColumn((*coliter).first);
      }
      catch (SQLException &e)
      {
        const string &exceptionmessage = e.getMessage();
        if (exceptionmessage.find("not found in result set") != string::npos)
        {
          cout << "Column " << (*coliter).first << " does not exist in "
               << tablename << ", creating it" << endl;
          cmd.str("");
          cmd << "ALTER TABLE "
              << tablename
              << " ADD "
              << (*coliter).first
              << " "
              << (*coliter).second;
          try
          {
            Statement *stmtup = DBconnection->createStatement();
            stmtup->executeUpdate(cmd.str());
          }
          catch (SQLException &e1)
          {
            cout << PHWHERE << " Exception caught: " << e1.getMessage() << endl;
          }
        }
      }
    }
    delete rs;
  }
  return 0;
}

int OnCalServer::add_calibrator_to_statustable(const std::string &calibratorname)
{
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  if (check_calibrator_in_statustable(calibratorname) == 0)
  {
    return 0;
  }
  const string &calibname = calibratorname;
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  cmd.str("");
  cmd << "ALTER TABLE " << successTable << " ADD COLUMN "
      << calibname << " int";
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    cout << "cmd: " << cmd.str() << endl;
    exit(1);
  }
  cmd.str("");
  cmd << "ALTER TABLE " << successTable << " ALTER COLUMN "
      << calibname << " SET DEFAULT " << OnCalDBCodes::INIT;
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    cout << "cmd: " << cmd.str() << endl;
    exit(1);
  }
  cmd.str("");
  cmd << "UPDATE " << successTable << " SET "
      << calibname << " = " << OnCalDBCodes::INIT;
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    cout << "cmd: " << cmd.str() << endl;
    exit(1);
  }

  return 0;
}

int OnCalServer::check_calibrator_in_statustable(const std::string &calibratorname)
{
  // replace this contraption by this sql command which returns 1 row if column exists
  // select * from information_schema.columns where table_name = 'oncal_status' and column_name = 'svxstripdeadmapcal';
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  string calibname = calibratorname;
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  cmd << "SELECT * FROM " << successTable << " LIMIT 1" << ends;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    cout << "Table " << successTable << " does not exist, your logic is off" << endl;
    exit(1);
  }
  ResultSetMetaData *meta = rs->getMetaData();
  unsigned int nocolumn = rs->getMetaData()->getColumnCount();
  // column names are lower case only, so convert string to lowercase
  // The bizarre cast here is needed for newer gccs
  transform(calibname.begin(), calibname.end(), calibname.begin(), (int (*)(int)) tolower);

  for (unsigned int i = 1; i <= nocolumn; i++)
  {
    if (meta->getColumnName(i) == calibname)
    {
      if (Verbosity() > 0)
      {
        cout << calibname << " is in " << successTable << endl;
      }
      return 0;
    }
  }
  // if we get here, the calibrator is not yet in the table
  delete rs;
  return -1;
}

int OnCalServer::check_create_successtable(const std::string &tablename)
{
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  cmd << "SELECT runnumber FROM " << tablename << " LIMIT 1" << ends;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << tablename << " does not exist, will create it" << endl;
    //      cout << "Message: " << e.getMessage() << endl;
  }
  if (!rs)
  {
    cmd.str("");
    cmd << "CREATE TABLE " << tablename << "(runnumber int NOT NULL, "
        << "startvaltime bigint, "
        << "endvaltime bigint, "
        << "beginrun timestamp(0) with time zone, "
        << "endrun timestamp(0) with time zone, "
        << "comment text, "
        << "primary key(runnumber))";
    cout << cmd.str() << endl;
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (SQLException &e)
    {
      cout << "Error, Message: " << e.getMessage() << endl;
      //      cout << "Message: " << e.getMessage() << endl;
    }
  }
  return 0;
}

void OnCalServer::recordDataBase(const bool bookkeep)
{
  recordDB = bookkeep;
  if (recordDB)
  {
    check_create_successtable(successTable);
  }
  return;
}

void OnCalServer::BeginTimeStamp(const PHTimeStamp &TimeStp)
{
  beginTimeStamp = TimeStp;
  cout << "OnCalServer::BeginTimeStamp: Setting BOR TimeStamp to " << beginTimeStamp << endl;
}

void OnCalServer::EndTimeStamp(const PHTimeStamp &TimeStp)
{
  endTimeStamp = TimeStp;
  cout << "OnCalServer::EndTimeStamp: Setting EOR TimeStamp to " << endTimeStamp << endl;
}

PHTimeStamp *
OnCalServer::GetLastGoodRunTS(OnCal *calibrator, const int irun)
{
  PHTimeStamp *ts = nullptr;
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return ts;
  }
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  ostringstream subsystable;
  subsystable << "oncal" << calibrator->Name();
  cmd << "SELECT runnumber FROM " << successTable << " where runnumber < "
      << irun << " and "
      << calibrator->Name() << " > 0 order by runnumber desc limit 1";
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << subsystable.str() << " does not exist" << endl;
    return ts;
  }
  if (rs->next())
  {
    RunToTime *rt = RunToTime::instance();
    int oldrun = rs->getInt("runnumber");
    ts = rt->getBeginTime(oldrun);
    cout << "Getting previous good run, current run: " << irun
         << ", previous good run: " << oldrun
         << " began ";
    ts->print();
    cout << endl;
  }
  else
  {
    cout << PHWHERE << " No previous good run found for run " << irun << endl;
  }
  delete rs;
  return ts;
}

int OnCalServer::SyncCalibTimeStampsToOnCal(const OnCal *calibrator, const int commit)
{
  vector<string> caltab;
  calibrator->GetPdbCalTables(caltab);
  vector<string>::const_iterator iter;
  for (iter = caltab.begin(); iter != caltab.end(); ++iter)
  {
    cout << "dealing with table: " << *iter << endl;
    SyncCalibTimeStampsToOnCal(calibrator, *iter, commit);
  }
  return 0;
}

int OnCalServer::SyncCalibTimeStampsToOnCal(const OnCal *calibrator, const std::string &table, const int commit)
{
  string name = calibrator->Name();
  odbc::Connection *con = nullptr;
  odbc::Connection *concalib = nullptr;
  ostringstream cmd;
  try
  {
    con = DriverManager::getConnection(database.c_str(), "phnxrc", "");
  }
  catch (SQLException &e)
  {
    cout << "Cannot connect to " << database.c_str() << endl;
    cout << e.getMessage() << endl;
    return -1;
  }
  try
  {
    concalib = DriverManager::getConnection("oncal", "phnxrc", "");
  }
  catch (SQLException &e)
  {
    cout << "Cannot connect to "
         << "oncal" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }
  Statement *stmt = nullptr;
  try
  {
    stmt = con->createStatement();
  }
  catch (SQLException &e)
  {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }

  PreparedStatement *stmt1 = nullptr;
  ResultSet *rs1 = nullptr;
  try
  {
    cmd.str("");
    cmd << "SELECT * from " << table << " where startvaltime = ?";
    stmt1 = concalib->prepareStatement(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }

  PreparedStatement *stmtupd = nullptr;
  try
  {
    cmd.str("");
    cmd << "update " << table << " set endvaltime = ? where startvaltime = ?";
    stmtupd = concalib->prepareStatement(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }

  cmd.str("");
  cmd << "select * from "
      << successTable
      << " where "
      << name
      << " > 0";
  //      << " > 0 and runnumber < 150000";
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    return -1;
  }
  while (rs->next())
  {
    int run = rs->getInt("runnumber");
    int startticks = rs->getLong("startvaltime");
    int endticks = rs->getLong("endvaltime");
    // int status = rs->getInt(name);
    //       cout << "run: " << run
    // 	   << ", status: " << status
    // 	   << ", startticks: " << startticks
    // 	   << ", endticks: " << endticks << endl;
    stmt1->setInt(1, startticks);
    try
    {
      rs1 = stmt1->executeQuery();
    }
    catch (SQLException &e)
    {
      cout << "Message: " << e.getMessage() << endl;
      return -1;
    }
    int ionce = 0;
    int isproblem = 0;
    int calibendval = 0;
    while (rs1->next())
    {
      calibendval = rs1->getInt("endvaltime");
      if (endticks != rs1->getInt("endvaltime"))
      {
        if (!isproblem)
        {
          cout << "endvaltime problem with run " << run << endl;
          cout << "endvaltime from oncal_status: " << endticks << endl;
          cout << "startvaltime from oncal_status: " << startticks << endl;
          cout << "endvaltime from calibrations DB: " << rs1->getInt("endvaltime") << endl;
          if (endticks < rs1->getInt("endvaltime"))
          {
            cout << "ENDTICKS smaller CALIB" << endl;
            //                      return -1;
          }
        }
        isproblem = 1;
      }
      else
      {
        if (isproblem)
        {
          cout << "endvaltime changes, check run " << run << endl;
          //                  return -1;
        }
      }
      // 	   cout << "starttime: " << rs1->getInt("startvaltime") << endl;
      // 	   cout << "endtime: " << rs1->getInt("endvaltime") << endl;
      ionce++;
    }
    if (isproblem)
    {
      cout << "Adjusting run " << run << endl;
      cout << "changing endvaltime from " << calibendval
           << " to " << endticks << endl;
      if (commit)
      {
        stmtupd->setInt(1, endticks);
        stmtupd->setInt(2, startticks);
        stmtupd->executeUpdate();
      }
    }
    if (!ionce)
    {
      cout << "Run " << run << " not found" << endl;
    }
    delete rs1;
  }
  delete rs;
  delete con;
  delete concalib;
  return 0;
}

int OnCalServer::SyncOncalTimeStampsToRunDB(const int commit)
{
  odbc::Connection *con = nullptr;
  RunToTime *rt = RunToTime::instance();
  ostringstream cmd;
  try
  {
    con = DriverManager::getConnection(database.c_str(), "phnxrc", "");
  }
  catch (SQLException &e)
  {
    cout << "Cannot connect to " << database.c_str() << endl;
    cout << e.getMessage() << endl;
    return -1;
  }
  Statement *stmt = nullptr;
  try
  {
    stmt = con->createStatement();
  }
  catch (SQLException &e)
  {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }

  PreparedStatement *stmtupd = nullptr;
  try
  {
    cmd.str("");
    cmd << "UPDATE oncal_status set endvaltime = ? where runnumber = ?";
    stmtupd = con->prepareStatement(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
  }

  cmd.str("");
  cmd << "select * from "
      << successTable;  //<< " where runnumber > 160000";
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Message: " << e.getMessage() << endl;
    return -1;
  }
  while (rs->next())
  {
    int run = rs->getInt("runnumber");
    int startticks = rs->getLong("startvaltime");
    int endticks = rs->getLong("endvaltime");
    int rtstartticks = 0;
    int rtendticks = 0;
    PHTimeStamp *rtstart = rt->getBeginTime(run);
    PHTimeStamp *rtend = rt->getEndTime(run);
    if (rtstart)
    {
      rtstartticks = rtstart->getTics();
      delete rtstart;
    }
    if (rtend)
    {
      rtendticks = rtend->getTics();
      delete rtend;
    }
    if (rtstartticks != startticks)
    {
      cout << "Run " << run
           << ": Start mismatch, oncal: " << startticks
           << ", rt: " << rtstartticks << endl;
    }
    if (rtendticks != endticks)
    {
      // exclude starttime=endtime in runtotime (some crashed calibrations can do this)
      // in this case the calibration adds 1 sec to starttime
      if (rtstartticks != rtendticks)
      {
        cout << "Run " << run
             << ": End mismatch, oncal: " << endticks
             << ", rt: " << rtendticks << endl;
        if (endticks > rtendticks)
        {
          cout << "BAD: endticks: " << endticks
               << ", rtendticks: " << rtendticks
               << endl;
          return -1;
        }
        if (commit)
        {
          stmtupd->setLong(1, rtendticks);
          stmtupd->setLong(2, run);
          stmtupd->executeUpdate();
        }
      }
      else
      {
        if (startticks != endticks - 1)
        {
          cout << "Run " << run
               << ": Start/End mismatch, Start: " << startticks
               << ", End: " << endticks << endl;
          endticks = startticks + 1;
          if (commit)
          {
            stmtupd->setLong(1, endticks);
            stmtupd->setLong(2, run);
            stmtupd->executeUpdate();
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            cout << "run " << run << " was twiddled by OnCal" << endl;
          }
        }
      }
    }
    //       cout << "run: " << run
    // 	   << ", status: " << status
    // 	   << ", startticks: " << startticks
    // 	   << ", endticks: " << endticks << endl;
  }
  delete rs;
  delete con;
  return 0;
}

int OnCalServer::CopyTables(const OnCal *calibrator, const int FromRun, const int ToRun, const int commit) const
{
  int iret = calibrator->CopyTables(FromRun, ToRun, commit);
  return iret;
}

int OnCalServer::CreateCalibration(OnCal *calibrator, const int myrunnumber, const string &what, const int commit)
{
  int iret = -1;
  runNum = myrunnumber;
  SetBorTime(myrunnumber);
  SetEorTime(myrunnumber);
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  add_calibrator_to_statustable(calibrator->Name());
  string table = "OnCal";
  table += calibrator->Name();
  check_create_subsystable(table);
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;

  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << myrunnumber;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  if (!rs->next())
  {
    insertRunNumInDB(successTable, myrunnumber);
  }
  delete rs;
  cmd.str("");
  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << myrunnumber << " and "
      << calibrator->Name() << " <= 0";
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << PHWHERE << " Exception caught, Message: "
         << e.getMessage() << endl;
    return -1;
  }
  if (rs->next() || testmode)
  {
    PdbBankManager *bankManager = PdbBankManager::instance();
    PdbApplication *application = bankManager->getApplication();
    application->setDBName("oncal");
    string tablecomment = "Subsytem provided";
    iret = calibrator->CreateCalibration(runnumber, what, tablecomment, commit);
    if (!iret)
    {
      cout << "Comment: " << tablecomment << endl;
      cout << "updating oncal status tables for " << runnumber << endl;
      if (commit)
      {
        CreateCalibrationUpdateStatus(calibrator, table, tablecomment, OnCalDBCodes::SUBSYSTEM);
      }
    }
    else
    {
      cout << "Calibratior " << calibrator->Name() << " for run " << runnumber << " failed" << endl;
      if (commit)
      {
        CreateCalibrationUpdateStatus(calibrator, table, tablecomment, OnCalDBCodes::FAILED);
      }
    }
  }
  else
  {
    cout << PHWHERE << " Run " << runnumber << " is already successfully calibrated for "
         << calibrator->Name() << endl;
  }
  return iret;
}

void OnCalServer::CreateCalibrationUpdateStatus(OnCal *calibrator, const string &table, const string &tablecomment, const int dbcode)
{
  updateDB(successTable.c_str(), calibrator->Name(), dbcode);
  insertRunNumInDB(table, RunNumber());
  updateDB(table, "comment", tablecomment, RunNumber(), true);
  ostringstream stringarg;
  stringarg.str("");
  stringarg << calibrator->CommitedToPdbCalOK();
  updateDB(table, "committed", stringarg.str(), RunNumber());
  stringarg.str("");
  stringarg << calibrator->VerificationOK();
  updateDB(table, "verified", stringarg.str(), RunNumber());
  odbc::Timestamp stp(time(nullptr));
  updateDB(table, "date", stp.toString(), RunNumber());
  time_t beginticks = beginTimeStamp.getTics();
  stringarg.str("");
  stringarg << beginticks;
  updateDB(table, "startvaltime", stringarg.str(), RunNumber());
  stp.setTime(beginticks);
  updateDB(table, "begintime", stp.toString(), RunNumber());
  time_t endticks = endTimeStamp.getTics();
  stringarg.str("");
  stringarg << endticks;
  updateDB(table, "endvaltime", stringarg.str(), RunNumber());
  stp.setTime(endticks);
  updateDB(table, "endtime", stp.toString(), RunNumber());
  updateDB(table, "cvstag", cvstag, RunNumber());
  vector<string> flist = calibrator->GetLocalFileList();
  if (flist.size())
  {
    string filelist = "";
    BOOST_FOREACH (string infile, flist)
    {
      filelist += infile;
      filelist += " ";
    }
    filelist.pop_back();  // strip empty space at end from loop
    cout << "FileList: " << filelist << endl;
    updateDB(table, "files", filelist, RunNumber());
  }
  return;
}

int OnCalServer::CopySnglTable(const string &pdbclass, const string &tablename, const int bankid, const int FromRun, const int ToRun, const int commit)
{
  int iret = CopySnglTableNewBankId(pdbclass, tablename, bankid, bankid, FromRun, ToRun, commit);
  return iret;
}

int OnCalServer::CopySnglTableNewBankId(const string &pdbclass, const string &tablename, const int bankid, const int Tobankid, const int FromRun, const int ToRun, const int commit)
{
  RunToTime *rt = RunToTime::instance();
  PHTimeStamp *ts = rt->getBeginTime(FromRun);
  assert(ts);
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  application->setDBName("oncal");
  PdbBankID BankID(bankid);
  PdbCalBank *pdbBank = bankManager->fetchBank(pdbclass.c_str(), BankID, tablename.c_str(), *ts);
  assert(pdbBank);
  PHTimeStamp *tstartnew = rt->getBeginTime(ToRun);
  PHTimeStamp *tendnew = rt->getEndTime(ToRun);
  application->startUpdate();
  PdbCalBank *pdbBanknew = (PdbCalBank *) (pdbBank->clone());
  ostringstream newdesc;
  newdesc << "copied from run " << FromRun;
  if (pdbBanknew)
  {
    pdbBanknew->setLength(pdbBank->getLength());
    if (Verbosity() > 0)
    {
      for (unsigned int i = 0; i < pdbBank->getLength(); i++)
      {
        cout << "orig: " << endl;
        pdbBank->printEntry(i);
        cout << "new: " << endl;
        pdbBanknew->printEntry(i);
      }
    }
    PHTimeStamp tsins;
    pdbBanknew->setInsertTime(tsins);
    pdbBanknew->setStartValTime(*tstartnew);
    pdbBanknew->setEndValTime(*tendnew);
    pdbBanknew->setDescription(newdesc.str().c_str());
    if (bankid != Tobankid)
    {
      pdbBanknew->setBankID(Tobankid);
    }
    if (Verbosity() > 0)
    {
      cout << "StartValTime: orig " << pdbBank->getStartValTime()
           << ", new: " << pdbBanknew->getStartValTime() << endl;
      cout << "EndValTime: orig " << pdbBank->getEndValTime()
           << ", new: " << pdbBanknew->getEndValTime() << endl;
    }
    if (commit)
    {
      application->commit(pdbBanknew);
    }
    else
    {
      application->abort();
    }
  }
  delete pdbBank;

  delete ts;
  delete tstartnew;
  delete tendnew;
  return 0;
}

int OnCalServer::ClosestGoodRun(OnCal *calibrator, const int irun, const int previous)
{
  RunToTime *rt = RunToTime::instance();
  PHTimeStamp *ts = rt->getBeginTime(irun);
  if (!ts)
  {
    cout << PHWHERE << "Unknown Run " << irun << endl;
    return -1;
  }
  int curstart = ts->getTics();
  delete ts;
  ts = rt->getEndTime(irun);
  int curend = curstart;
  if (ts)
  {
    curend = ts->getTics();
    delete ts;
  }
  int closestrun = -1;
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;

  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber < "
      << irun << " and "
      << calibrator->Name() << " = 1 order by runnumber desc limit 1";
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  int prevrun = -1;
  unsigned int prevend = 0;
  if (rs->next())
  {
    prevrun = rs->getInt("runnumber");
    unsigned int prevstart = rs->getLong("startvaltime");
    prevend = rs->getLong("endvaltime");
    cout << "previous run: " << prevrun
         << ", start: " << prevstart
         << ", end: " << prevend
         << endl;
  }
  else
  {
    if (Verbosity() > 0)
    {
      cout << PHWHERE << " No previous good run found for run " << irun << endl;
    }
  }
  delete rs;
  closestrun = prevrun;
  if (previous == fetchrun::PREVIOUS)
  {
    if (Verbosity() > 0)
    {
      cout << "Closest previous run is " << closestrun << endl;
    }
    return closestrun;
  }
  cmd.str("");
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber > "
      << irun << " and "
      << calibrator->Name() << " = 1 order by runnumber asc limit 1";
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  int nextrun = -1;
  unsigned int nextstart = 0;
  if (rs->next())
  {
    nextrun = rs->getInt("runnumber");
    nextstart = rs->getLong("startvaltime");
    unsigned int nextend = rs->getLong("endvaltime");
    if (Verbosity() > 0)
    {
      cout << "next run: " << nextrun
           << ", start: " << nextstart
           << ", end: " << nextend
           << endl;
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      cout << PHWHERE << " No next good run found for run " << irun << endl;
    }
  }
  delete rs;
  int tdiffprev = curstart - prevend;
  int tdiffnext;
  if (nextstart > 0)
  {
    tdiffnext = nextstart - curend;
  }
  else
  {
    // just make it larger then previous run time diff
    tdiffnext = tdiffprev + 1;
  }
  if (Verbosity() > 0)
  {
    cout << "diff prev: " << tdiffprev
         << ", next: " << tdiffnext
         << endl;
  }
  if (tdiffprev < tdiffnext)
  {
    closestrun = prevrun;
  }
  else
  {
    closestrun = nextrun;
  }
  if (Verbosity() > 0)
  {
    cout << "closest run: " << closestrun << endl;
  }
  return closestrun;
}

int OnCalServer::OverwriteCalibration(OnCal *calibrator, const int runno, const int commit, const int FromRun)
{
  if (FromRun < 0)
  {
    return -1;
  }
  int iret = CopyTables(calibrator, FromRun, runno, commit);
  return iret;
}

int OnCalServer::FixMissingCalibration(OnCal *calibrator, const int runno, const int commit, const int fromrun)
{
  int iret = -1;
  // find this run in oncal_status
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  runNum = runno;
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;

  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << runno;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  if (!rs->next())
  {
    insertRunNumInDB(successTable, runNum);
  }
  delete rs;
  cmd.str("");
  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << runno << " and "
      << calibrator->Name() << " <= 0";
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << PHWHERE << " Exception caught, Message: "
         << e.getMessage() << endl;
    return -1;
  }
  if (rs->next())
  {
    int FromRun;
    if (fromrun > 0)
    {
      FromRun = fromrun;
    }
    else
    {
      FromRun = ClosestGoodRun(calibrator, runno);
      if (FromRun < 0)
      {
        cout << "ClosestGoodRun returned bad runnumber: " << FromRun << endl;
        return -1;
      }
    }
    cout << "Going to copy calibration for run " << runno
         << " from run " << FromRun << endl;

    iret = OverwriteCalibration(calibrator, runno, commit, FromRun);
    if (!iret)
    {
      int newstatus = 0;
      if (FromRun < runno)
      {
        newstatus = OnCalDBCodes::COPIEDPREVIOUS;
      }
      else
      {
        newstatus = OnCalDBCodes::COPIEDLATER;
      }
      string table = "OnCal";
      table += calibrator->Name();
      ostringstream comment;
      comment << " CopiedRun(" << FromRun << ")";
      cout << "updating oncal status tables for " << runno << endl;
      if (commit)
      {
        updateDB(successTable.c_str(), calibrator->Name(), newstatus);
        insertRunNumInDB(table, runNum);
        updateDB(table, "comment", comment.str(), runNum, true);
        updateDB(table.c_str(), "committed", true);
      }
    }
  }
  else
  {
    cout << "Run " << runno
         << " has a good calibrations, doing nothing" << endl;
  }
  delete rs;
  return iret;
}

int OnCalServer::SetBorTime(const int runno)
{
  //  recoConsts *rc = recoConsts::instance();
  RunToTime *runTime = RunToTime::instance();

  PHTimeStamp *BorTimeStp(runTime->getBeginTime(runno));
  if (!BorTimeStp)
  {
    cout << PHWHERE << "Cannot get begin time for run " << runno << endl;
    cout << "Exiting" << endl;
    exit(1);
  }
  BeginTimeStamp(*BorTimeStp);

  // enter begin run timestamp into rc flags
  PHTimeStamp BeginRunTimeStamp(*BorTimeStp);
  //  rc->set_TimeStamp(BeginRunTimeStamp);
  cout << "OnCalServer::SetBorTime from RunToTime was found for run : " << runno << " to ";
  BeginRunTimeStamp.print();
  cout << endl;

  delete BorTimeStp;
  return 0;
}

int OnCalServer::SetEorTime(const int runno)
{
  //  recoConsts *rc = recoConsts::instance();
  RunToTime *runTime = RunToTime::instance();

  time_t eorticks = 0;

  time_t borticks = 0;  //(rc->get_TimeStamp()).getTics();
  PHTimeStamp *EorTimeStp(runTime->getEndTime(runno));
  if (EorTimeStp)
  {
    eorticks = EorTimeStp->getTics();
  }
  else
  {
    EorTimeStp = new PHTimeStamp(eorticks);
  }
  // if end of run timestamp missing or smaller-equal borstamp eor = bor+1 sec
  if (eorticks <= borticks)
  {
    eorticks = borticks + 1;
    EorTimeStp->setTics(eorticks);
  }
  EndTimeStamp(*EorTimeStp);
  cout << "OnCalServer::SetEorTime: setting eor time to ";
  EorTimeStp->print();
  cout << endl;
  delete EorTimeStp;
  return 0;
}

int OnCalServer::GetRunTimeTicks(const int runno, time_t &borticks, time_t &eorticks)
{
  RunToTime *runTime = RunToTime::instance();
  PHTimeStamp *TimeStp(runTime->getBeginTime(runno));
  if (!TimeStp)
  {
    cout << PHWHERE << "Cannot get begin time for run " << runno << endl;
    cout << "Exiting" << endl;
    exit(1);
  }
  borticks = TimeStp->getTics();
  delete TimeStp;
  TimeStp = runTime->getEndTime(runno);
  if (TimeStp)
  {
    eorticks = TimeStp->getTics();
    delete TimeStp;
  }
  else
  {
    eorticks = 0;
  }
  // if end of run timestamp missing or smaller-equal borstamp eor = bor+1 sec
  if (eorticks <= borticks)
  {
    eorticks = borticks + 1;
  }
  return 0;
}

int OnCalServer::requiredCalibration(SubsysReco *reco, const std::string &calibratorname)
{
  map<string, set<SubsysReco *> >::iterator iter;
  if (check_calibrator_in_statustable(calibratorname))
  {
    cout << PHWHERE << " the calibrator " << calibratorname << " is unknown to me" << endl;
    return -1;
  }
  iter = requiredCalibrators.find(calibratorname);
  if (iter != requiredCalibrators.end())
  {
    iter->second.insert(reco);
  }
  else
  {
    set<SubsysReco *> subsys;
    subsys.insert(reco);
    requiredCalibrators[calibratorname] = subsys;
  }
  return 0;
}

int OnCalServer::FindClosestCalibratedRun(const int irun)
{
  RunToTime *rt = RunToTime::instance();
  PHTimeStamp *ts = rt->getBeginTime(irun);
  if (!ts)
  {
    cout << PHWHERE << "Unknown Run " << irun << endl;
    return -1;
  }
  if (requiredCalibrators.size() == 0)
  {
    cout << PHWHERE << "No required calibrations given" << endl;
    return irun;
  }
  int curstart = ts->getTics();
  delete ts;
  ts = rt->getEndTime(irun);
  int curend = curstart;
  if (ts)
  {
    curend = ts->getTics();
    delete ts;
  }
  int closestrun = -1;
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -1;
  }
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;
  map<string, set<SubsysReco *> >::const_iterator iter;
  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber <= "
      << irun;
  for (iter = requiredCalibrators.begin(); iter != requiredCalibrators.end(); ++iter)
  {
    cmd << " and " << iter->first << " > 0 ";
  }

  cmd << " order by runnumber desc limit 1";
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  int prevrun = 0;
  unsigned int prevend = 0;
  if (rs->next())
  {
    prevrun = rs->getInt("runnumber");
    unsigned int prevstart = rs->getLong("startvaltime");
    prevend = rs->getLong("endvaltime");
    if (prevrun != irun)
    {
      cout << "previous run: " << prevrun
           << ", start: " << prevstart
           << ", end: " << prevend
           << endl;
    }
  }
  else
  {
    cout << PHWHERE << " No previous good run found for run " << irun << endl;
  }
  delete rs;
  // if the current run fullfills requirements return immediately
  if (prevrun == irun)
  {
    cout << "closest run with required calibs is current run: " << irun << endl;
    return irun;
  }
  cmd.str("");
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber > "
      << irun;
  for (iter = requiredCalibrators.begin(); iter != requiredCalibrators.end(); ++iter)
  {
    cmd << " and " << iter->first << " > 0 ";
  }

  cmd << " order by runnumber asc limit 1";
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -1;
  }
  int nextrun = 0;
  unsigned int nextstart = 0;
  if (rs->next())
  {
    nextrun = rs->getInt("runnumber");
    nextstart = rs->getLong("startvaltime");
    unsigned int nextend = rs->getLong("endvaltime");
    cout << "next run: " << nextrun
         << ", start: " << nextstart
         << ", end: " << nextend
         << endl;
  }
  else
  {
    cout << PHWHERE << " No next good run found for run " << irun << endl;
  }
  delete rs;
  int tdiffprev = curstart - prevend;
  int tdiffnext;
  if (nextstart > 0)
  {
    tdiffnext = nextstart - curend;
  }
  else
  {
    // just make it larger then previous run time diff
    tdiffnext = tdiffprev + 1;
  }
  if (tdiffprev < tdiffnext)
  {
    closestrun = prevrun;
  }
  else
  {
    closestrun = nextrun;
  }
  cout << "closest run with required calibs: " << closestrun << endl;
  return closestrun;
}

int OnCalServer::FillRunListFromFileList()
{
  // get sync managers and ask them for input managers and ask those for
  // their respective file lists. BOOST_FOREACH allows us to do loops without
  // having to explicietly declare all those vectors/lists
  BOOST_FOREACH (Fun4AllSyncManager *sync, SyncManagers)
  {
    BOOST_FOREACH (Fun4AllInputManager *inmgr, sync->GetInputManagers())
    {
      BOOST_FOREACH (string infile, inmgr->GetFileList())
      {
        pair<int, int> runseg = Fun4AllUtils::GetRunSegment(infile);
        runlist.insert(runseg.first);
      }
    }
  }
  return 0;
}

int OnCalServer::AdjustRichTimeStampForMultipleRuns()
{
  int firstrun = *runlist.begin();
  int lastrun = *runlist.rbegin();
  time_t dummy;
  time_t beginticks;
  time_t endticks;
  string table = "OnCalRichCal";
  check_create_subsystable(table);
  GetRunTimeTicks(firstrun, beginticks, dummy);
  GetRunTimeTicks(lastrun, dummy, endticks);
  ostringstream stringarg;
  stringarg << OnCalDBCodes::COVERED;
  //  set<int>::const_iterator runiter;
  /*
    for (runiter = runlist.begin(); runiter != runlist.end(); runiter++)
    {
    updateDB(successTable, "RichCal", stringarg.str(), *runiter);
    }
    stringarg.str("");
    stringarg << OnCalDBCodes::SUCCESS;

    updateDB(successTable, "RichCal", stringarg.str(), firstrun);
  */
  odbc::Timestamp stp;
  stringarg.str("");
  stringarg << beginticks;
  updateDB(table, "startvaltime", stringarg.str(), firstrun);
  stp.setTime(beginticks);
  updateDB(table, "begintime", stp.toString(), firstrun);
  stringarg.str("");
  stringarg << endticks;
  updateDB(table, "endvaltime", stringarg.str(), firstrun);
  stp.setTime(endticks);
  updateDB(table, "endtime", stp.toString(), firstrun);
  /*
    string tablename = "calibrichadc";
    odbc::Connection *con = 0;
    ostringstream cmd;
    try
    {
    con = DriverManager::getConnection("oncal", "phnxrc", "");
    }
    catch (SQLException& e)
    {
    cout << "Cannot connect to " << database.c_str() << endl;
    cout << e.getMessage() << endl;
    return -1;
    }
    Statement *stmt = 0;
    Statement *stmtup = 0;
    try
    {
    stmt = con->createStatement();
    stmtup = con->createStatement();
    }
    catch (SQLException& e)
    {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
    }

    ResultSet *rs1 = 0;
    cmd.str("");
    cmd << "SELECT  endvaltime from " << tablename
    << " where bankid = 1 and startvaltime = " << beginticks;
    cout << "sql cmd: " << cmd.str() << endl;
    try
    {
    rs1 = stmt->executeQuery(cmd.str());
    }
    catch (SQLException& e)
    {
    cout << "Cannot create statement" << endl;
    cout << e.getMessage() << endl;
    return -1;
    }
    if (rs1->next())
    {
    cout << "Endcaltime: " << rs1->getInt("endvaltime") << endl;
    cout << "future endvaltime: " << endticks << endl;
    cmd.str("");
    cmd << "Update " << tablename
    << " set endvaltime = " << endticks
    << " where bankid = 1 and startvaltime = "
    << beginticks;
    stmtup->executeUpdate(cmd.str());

    }
    else
    {
    cout << "Could not find startvaltime " << beginticks
    << "from run " << firstrun << endl;
    }

  */

  return 0;
}

int OnCalServer::GetCalibStatus(const string &calibname, const int runno)
{
  int iret = -3;
  if (!connectDB())
  {
    cout << "could not connect to " << database << endl;
    return -4;
  }
  Statement *stmt = DBconnection->createStatement();
  ostringstream cmd;

  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT " << calibname << " FROM "
      << successTable << " where runnumber = "
      << runno;
  cout << "exec " << cmd.str() << endl;
  ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (SQLException &e)
  {
    cout << "Table " << successTable << " does not exist" << endl;
    return -5;
  }
  if (rs->next())
  {
    iret = rs->getInt(calibname.c_str());
  }
  else
  {
    cout << PHWHERE << " No calib status for " << calibname
         << " for " << runno << endl;
  }
  delete rs;
  return iret;
}

void OnCalServer::TestMode(const int i)
{
  const char *logname = getenv("LOGNAME");
  if (logname)
  {
    if (strcmp(logname, "phnxcal") == 0 || strcmp(logname, "anatrain") == 0)
    {
      cout << "phnxcal,anatrain account is not allowed to run in testmode" << endl;
    }
    else
    {
      testmode = i;
    }
  }
  else
  {
    cout << "could not get account via env var LOGNAME, not setting testmode" << endl;
  }
  return;
}

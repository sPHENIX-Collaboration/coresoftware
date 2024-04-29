#include "OnCalServer.h"
#include "OnCal.h"
#include "OnCalDBCodes.h"
#include "OnCalHistoBinDefs.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>  // for Fun4AllServer, Fun4AllServe...
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHTimeStamp.h>  // for PHTimeStamp, operator<<
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/RunToTime.h>

#include <RtypesCore.h>  // for Stat_t
#include <TDirectory.h>  // for TDirectoryAtomicAdapter
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>  // for TNamed
#include <TROOT.h>
#include <TString.h>  // for TString

// odbc++ classes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/preparedstatement.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException, Timestamp
#pragma GCC diagnostic pop

#include <unistd.h>
#include <algorithm>
#include <cassert>  // for assert
#include <cctype>   // for tolower
#include <cstdlib>
#include <cstring>  // for strcmp
#include <fstream>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <list>
#include <map>
#include <sstream>
#include <utility>  // for pair

const static std::string cvstag = "OnCalv86";

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
    std::cout << PHWHERE << "Screwup - the end validity time is not set" << std::endl;
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
    std::cout << PHWHERE << "Screwup - the begin validity time is not set" << std::endl;
    exit(1);
  }
}
//---------------------------------------------------------------------

void OnCalServer::dumpHistos()
{
  std::ostringstream filename;
  std::string fileprefix = "./";

  if (getenv("ONCAL_SAVEDIR"))
  {
    fileprefix = getenv("ONCAL_SAVEDIR");
    fileprefix += "/";
  }

  int compress = 3;
  std::map<std::string, std::set<std::string> >::const_iterator iter;
  //  std::map<std::string, TH1 *>::const_iterator hiter;
  TH1 *histo;
  std::set<std::string>::const_iterator siter;
  for (iter = calibratorhistomap.begin(); iter != calibratorhistomap.end(); ++iter)
  {
    filename.str("");
    filename << fileprefix << "Run_"
             << RunNumber()
             << "_" << iter->first << ".root";
    TFile *hfile = new TFile(filename.str().c_str(), "RECREATE",
                             "Created by Online Calibrator", compress);
    std::cout << "OnCalServer::dumpHistos() Output root file: " << filename.str() << std::endl;
    for (siter = (iter->second).begin(); siter != (iter->second).end(); ++siter)
    {
      histo = dynamic_cast<TH1 *>(getHisto(*siter));
      if (histo)
      {
        histo->Write();
      }
      else
      {
        std::cout << PHWHERE << "Histogram "
                  << *siter << " not found, will not be saved in "
                  << filename.str() << std::endl;
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
    std::string calibratorname = Calibrator->Name();
    std::map<std::string, std::set<std::string> >::iterator iter;
    iter = calibratorhistomap.find(calibratorname);
    if (iter != calibratorhistomap.end())
    {
      (iter->second).insert(h1d->GetName());
    }
    else
    {
      std::set<std::string> newset;
      newset.insert(h1d->GetName());
      newset.insert("OnCalServerVars");
      calibratorhistomap[calibratorname] = newset;
    }
  }
  Fun4AllServer::registerHisto(h1d, replace);
  return;
}

void OnCalServer::unregisterHisto(const std::string &calibratorname)
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
    std::cout << nEvents << " events, testing" << std::endl;
    unsigned int j = 0;
    unsigned int ical = 0;
    std::vector<std::pair<SubsysReco *, PHCompositeNode *> >::const_iterator iter;
    for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
    {
      OnCal *oncal = dynamic_cast<OnCal *>(iter->first);
      if (oncal)
      {
        ical++;
        std::cout << "Name: " << oncal->Name()
                  << " is " << oncal->AllDone() << std::endl;
        j += oncal->AllDone();
      }
    }
    if (j == ical)
    {
      std::cout << "Everyone is done after "
                << nEvents << " Events" << std::endl;
      i = 1;
    }
  }
  return i;
}

int OnCalServer::BeginRun(const int runno)
{
  if (runno <= 0)
  {
    std::cout << PHWHERE << "Invalid Run Number: " << runno << std::endl;
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
  std::vector<std::pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  // copy the subsys reco pointers to another set for
  // easier search (we only need the pointers to find
  // the subsystems with special timestamp/runnumber needs
  std::set<SubsysReco *> NeedOtherTimeStamp;
  std::map<std::string, std::set<SubsysReco *> >::const_iterator miter;
  std::set<SubsysReco *>::const_iterator siter;
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
  std::string currdir = gDirectory->GetPath();
  std::set<std::string> droplist;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    std::ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << (*iter).second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }
    OnCal *oncal = dynamic_cast<OnCal *>((*iter).first);
    if (oncal)
    {
      std::string table = "OnCal";
      table += (*iter).first->Name();
      check_create_subsystable(table);
      insertRunNumInDB(table, runNum);
      std::string calibname = (*iter).first->Name();
      add_calibrator_to_statustable(calibname);
      std::set<int>::const_iterator runiter;
      int calibstatus = GetCalibStatus(calibname, runNum);
      if (calibstatus > 0 && testmode == 0)
      {
        std::cout << calibname << " already ran for run " << runNum << std::endl;
        droplist.insert(calibname);
        unregisterSubsystem(oncal);
        unregisterHisto(calibname);
      }
      else
      {
        std::ostringstream stringarg;
        stringarg << OnCalDBCodes::STARTED;
        for (runiter = runlist.begin(); runiter != runlist.end(); ++runiter)
        {
          updateDB(successTable, calibname, stringarg.str(), *runiter);
        }
      }
    }
    if (NeedOtherTimeStamp.find((*iter).first) != NeedOtherTimeStamp.end())
    {
      std::cout << "changing timestamp for " << (*iter).first->Name() << std::endl;
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
        std::cout << PHWHERE << "Module " << (*iter).first->Name() << " issued Abort Run, exiting" << std::endl;
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
    std::cout << "No Events read, you probably gave me an empty filelist" << std::endl;
    return -1;
  }
  PdbApplication *app = PdbApplication::instance();
  app->setDBName("oncal");
  int i = 0;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();

  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    std::ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << (*iter).second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }
    else
    {
      if (Verbosity() > 2)
      {
        std::cout << "End: cded to " << newdirname.str().c_str() << std::endl;
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
    std::ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << (*iter).second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }

    std::string CalibratorName = oncal->Name();

    int verificationstatus = oncal->VerificationOK();
    int databasecommitstatus = oncal->CommitedToPdbCalOK();

    // report success database the status of the calibration
    if (recordDB)
    {
      std::string table = "OnCal";
      table += CalibratorName;

      std::ostringstream stringarg;
      if (databasecommitstatus == OnCalDBCodes::SUCCESS)
      {
        stringarg << OnCalDBCodes::COVERED;
      }
      else
      {
        stringarg << OnCalDBCodes::FAILED;
      }
      std::set<int>::const_iterator runiter;
      for (runiter = runlist.begin(); runiter != runlist.end(); ++runiter)
      {
        updateDB(successTable, CalibratorName, stringarg.str(), *runiter);
      }
      // update the first run which was used in the calibration
      // with the real status
      updateDB(successTable, CalibratorName, databasecommitstatus);

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

      std::string filelist = "";
      for (Fun4AllSyncManager *sync : SyncManagers)
      {
        for (Fun4AllInputManager *inmgr : sync->GetInputManagers())
        {
          for (const std::string &infile : inmgr->GetFileOpenedList())
          {
            filelist += (infile).substr(((infile).find_last_of('/') + 1), (infile).size());
            filelist += " ";  // this needs to be stripped again for last entry
          }
        }
      }
      filelist.pop_back();  // strip empty space at end from loop
      std::cout << "FileList: " << filelist << std::endl;
      updateDB(table, "files", filelist, RunNumber());
      updateDB(table, "cvstag", cvstag, RunNumber());
    }

    std::cout << "SERVER SUMMARY: " << oncal->Name() << "  "
              << (verificationstatus == 1 ? "Verification: SUCCESS  " : "Verification: FAILURE  ")
              << (databasecommitstatus == 1 ? "DB commit: SUCCESS  " : "DB commit: FAILURE  ")
              << std::endl;

    printStamps();
  }
  gROOT->cd(currdir.c_str());
  dumpHistos();  // save the histograms in files
  return i;
}
//---------------------------------------------------------------------

void OnCalServer::Print(const std::string &what) const
{
  Fun4AllServer::Print(what);
  if (what == "ALL" || what == "CALIBRATOR")
  {
    // loop over the map and print out the content
    // (name and location in memory)

    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Calibrators in OnCalServer:" << std::endl;

    std::vector<std::pair<SubsysReco *, PHCompositeNode *> >::const_iterator miter;
    for (miter = Subsystems.begin();
         miter != Subsystems.end(); ++miter)
    {
      OnCal *oncal = dynamic_cast<OnCal *>((*miter).first);
      if (oncal)
      {
        std::cout << oncal->Name() << std::endl;
      }
    }
    std::cout << std::endl;
  }
  if (what == "ALL" || what == "REQUIRED")
  {
    // loop over the map and print out the content
    // (name and location in memory)

    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of required Calibrations in OnCalServer:" << std::endl;

    std::map<std::string, std::set<SubsysReco *> >::const_iterator iter;
    std::set<SubsysReco *>::const_iterator siter;
    for (iter = requiredCalibrators.begin();
         iter != requiredCalibrators.end(); ++iter)
    {
      std::cout << iter->first << " calibrations are needed by " << std::endl;
      for (siter = iter->second.begin(); siter != iter->second.end(); ++siter)
      {
        std::cout << (*siter)->Name() << std::endl;
      }
    }
    std::cout << std::endl;
  }
  if (what == "ALL" || what == "FILES")
  {
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of PRDF Files in OnCalServer:" << std::endl;
    for (Fun4AllSyncManager *sync : SyncManagers)
    {
      for (Fun4AllInputManager *inmgr : sync->GetInputManagers())
      {
        for (const std::string &infile : inmgr->GetFileList())
        {
          std::cout << "File: " << infile << std::endl;
        }
      }
    }
  }
  if (what == "ALL" || what == "RUNS")
  {
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Run Numbers in OnCalServer:" << std::endl;
    std::set<int>::const_iterator liter;
    for (liter = runlist.begin(); liter != runlist.end(); ++liter)
    {
      std::cout << "Run : " << *liter << std::endl;
    }
  }

  return;
}

void OnCalServer::printStamps()
{
  std::cout << std::endl
            << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << "*    VALIDITY RANGE FOR THIS CALIBRATION  *" << std::endl;
  std::cout << "*                                         *" << std::endl;
  std::cout << "* Used Run     :   ";
  std::cout << runNum << std::endl;
  std::cout << std::endl;
  std::cout << "* Begin Valid  :   ";
  beginTimeStamp.print();
  std::cout << std::endl;
  std::cout << "* End Valid    :   ";
  endTimeStamp.print();
  std::cout << std::endl;
  std::cout << "*                                         *" << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << std::endl
            << std::endl
            << std::endl;
}

//---------------------------------------------------------------------

void OnCalServer::RunNumber(const int runnum)
{
  runNum = runnum;
  SetBorTime(runnum);
  if (recordDB)
  {
    std::set<int>::const_iterator runiter;
    time_t beginrunticks;
    time_t endrunticks;
    std::ostringstream stringarg;
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
          odbc::DriverManager::getConnection(database.c_str(), "phnxrc", "");
    }
    catch (odbc::SQLException &e)
    {
      std::cout << "Cannot connect to " << database.c_str() << std::endl;
      std::cout << e.getMessage() << std::endl;
      std::cout << "countdown: " << countdown << std::endl;
      countdown--;
      failure = true;
      sleep(100);  // try again in 100 secs
    }
  }
  if (failure)
  {
    std::cout << "could not connect to DB after 10 tries in 1000 secs, giving up" << std::endl;
    exit(-1);
  }
  std::cout << "connected to " << database.c_str() << " database." << std::endl;
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

bool OnCalServer::insertRunNumInDB(const std::string &DBtable, const int runno)
{
  if (findRunNumInDB(DBtable, runno))
  {
    return true;
  }

  std::cout << "new row will be created in DB for run " << runno << std::endl;

  odbc::Statement *statement = nullptr;
  statement = DBconnection->createStatement();
  std::ostringstream cmd;
  cmd << "INSERT INTO "
      << DBtable
      << " (runnumber) VALUES ("
      << runno << ")";

  if (Verbosity() == 1)
  {
    std::cout << "in function OnCalServer::insertRunNumInDB() ... ";
    std::cout << "executing SQL statements ..." << std::endl;
    std::cout << cmd.str() << std::endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

bool OnCalServer::findRunNumInDB(const std::string &DBtable, const int runno)
{
  if (!DBconnection)
  {
    connectDB();
  }
  odbc::Statement *statement = nullptr;
  odbc::ResultSet *rs = nullptr;
  std::ostringstream cmd;
  cmd << "SELECT runnumber FROM "
      << DBtable
      << " WHERE runnumber = "
      << runno;

  statement = DBconnection->createStatement();

  if (Verbosity() == 1)
  {
    std::cout << "in function OnCalServer::findRunNumInDB() ";
    std::cout << "executing SQL statement ..." << std::endl
              << cmd.str() << std::endl;
  }

  try
  {
    rs = statement->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE << " exception caught: " << e.getMessage() << std::endl;
    return false;
  }

  int entry = 0;
  if (rs->next())
  {
    try
    {
      entry = rs->getInt("runnumber");
    }
    catch (odbc::SQLException &e)
    {
      std::cout << PHWHERE << " exception caught: " << e.getMessage() << std::endl;
      return false;
    }
  }
  else
  {
    return false;
  }
  std::cout << "run number " << entry << " already exists in DB" << std::endl;
  return true;
}

bool OnCalServer::updateDBRunRange(const std::string &table, const std::string &column, const int entry, const int firstrun, const int lastrun)
{
  if (!DBconnection)
  {
    connectDB();
  }

  odbc::Statement *statement = nullptr;

  std::string command = "UPDATE ";
  command += table;
  command += " SET ";
  command += column;
  command += " = ";
  command += std::to_string(entry);
  command += " WHERE runnumber >= ";
  command += std::to_string(firstrun);
  command += " and runnumber <= ";
  command += std::to_string(lastrun);

  if (Verbosity() == 1)
  {
    std::cout << "in function OnCalServer::updateDB() ... ";
    std::cout << "executin SQL statement ... " << std::endl;
    std::cout << command << std::endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command);
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

bool OnCalServer::updateDB(const std::string &table, const std::string &column, int entry)
{
  if (!DBconnection)
  {
    connectDB();
  }

  odbc::Statement *statement = nullptr;

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
    std::cout << "in function OnCalServer::updateDB() ... ";
    std::cout << "executin SQL statement ... " << std::endl;
    std::cout << command.Data() << std::endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command.Data());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
    return false;
  }

  return true;
}
//---------------------------------------------------------------------

bool OnCalServer::updateDB(const std::string &table, const std::string &column, bool entry)
{
  if (!DBconnection)
  {
    connectDB();
  }
  odbc::Statement *statement = nullptr;

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
    std::cout << "in function OnCalServer::updateDB() ... ";
    std::cout << "executin SQL statement ... " << std::endl;
    std::cout << command.Data() << std::endl;
  }
  statement = DBconnection->createStatement();

  try
  {
    statement->executeUpdate(command.Data());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
    return false;
  }

  return true;
}

//---------------------------------------------------------------------

int OnCalServer::updateDB(const std::string &table, const std::string &column,
                          const time_t ticks)
{
  if (!DBconnection)
  {
    connectDB();
  }
  odbc::Statement *statement = nullptr;

  std::ostringstream cmd;
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
    std::cout << "in function OnCalServer::updateDB() ... ";
    std::cout << "executin SQL statement ... " << std::endl;
    std::cout << cmd.str() << std::endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
    return -1;
  }
  return 0;
}
//---------------------------------------------------------------------

bool OnCalServer::updateDB(const std::string &table, const std::string &column,
                           const std::string &entry, const int runno, const bool append)
{
  if (!DBconnection)
  {
    connectDB();
  }

  odbc::Statement *statement = nullptr;

  statement = DBconnection->createStatement();

  std::string comment = "";
  std::ostringstream cmd;
  if (append)
  {
    odbc::ResultSet *rs = nullptr;
    std::ostringstream query;
    query << "SELECT * FROM "
          << table
          << " WHERE runnumber = "
          << runno;

    try
    {
      rs = statement->executeQuery(query.str());
    }
    catch (odbc::SQLException &e)
    {
      std::cout << "in function OnCalServer::updateDB() ... ";
      std::cout << "run number " << runno << "not found in DB" << std::endl;
      std::cout << e.getMessage() << std::endl;
    }

    rs->next();
    try
    {
      comment = rs->getString(column);
      comment += " ";  // add empty space between comments
    }
    catch (odbc::SQLException &e)
    {
      std::cout << "in function OnCalServer::updateDB() ... " << std::endl;
      std::cout << "nothing to append." << std::endl;
      std::cout << e.getMessage() << std::endl;
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
    std::cout << "in function OnCalServer::updateDB() ... ";
    std::cout << "executin SQL statement ... " << std::endl;
    std::cout << cmd.str() << std::endl;
  }

  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << e.getMessage() << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  std::vector<std::pair<std::string, std::string> > calibrator_columns;
  std::vector<std::pair<std::string, std::string> >::const_iterator coliter;
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

  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT * FROM " << tablename << " LIMIT 1" << std::ends;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << tablename << " does not exist, will create it" << std::endl;
    //      std::cout << "Message: " << e.getMessage() << std::endl;
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
      catch (odbc::SQLException &e)
      {
        const std::string &exceptionmessage = e.getMessage();
        if (exceptionmessage.find("not found in result set") != std::string::npos)
        {
          std::cout << "Column " << (*coliter).first << " does not exist in "
                    << tablename << ", creating it" << std::endl;
          cmd.str("");
          cmd << "ALTER TABLE "
              << tablename
              << " ADD "
              << (*coliter).first
              << " "
              << (*coliter).second;
          try
          {
            odbc::Statement *stmtup = DBconnection->createStatement();
            stmtup->executeUpdate(cmd.str());
          }
          catch (odbc::SQLException &e1)
          {
            std::cout << PHWHERE << " Exception caught: " << e1.getMessage() << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  if (check_calibrator_in_statustable(calibratorname) == 0)
  {
    return 0;
  }
  const std::string &calibname = calibratorname;
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  cmd.str("");
  cmd << "ALTER TABLE " << successTable << " ADD COLUMN "
      << calibname << " int";
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
    std::cout << "cmd: " << cmd.str() << std::endl;
    exit(1);
  }
  cmd.str("");
  cmd << "ALTER TABLE " << successTable << " ALTER COLUMN "
      << calibname << " SET DEFAULT " << OnCalDBCodes::INIT;
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
    std::cout << "cmd: " << cmd.str() << std::endl;
    exit(1);
  }
  cmd.str("");
  cmd << "UPDATE " << successTable << " SET "
      << calibname << " = " << OnCalDBCodes::INIT;
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
    std::cout << "cmd: " << cmd.str() << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  std::string calibname = calibratorname;
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT * FROM " << successTable << " LIMIT 1" << std::ends;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
    std::cout << "Table " << successTable << " does not exist, your logic is off" << std::endl;
    exit(1);
  }
  odbc::ResultSetMetaData *meta = rs->getMetaData();
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
        std::cout << calibname << " is in " << successTable << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT runnumber FROM " << tablename << " LIMIT 1" << std::ends;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << tablename << " does not exist, will create it" << std::endl;
    //      std::cout << "Message: " << e.getMessage() << std::endl;
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
    std::cout << cmd.str() << std::endl;
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException &e)
    {
      std::cout << "Error, Message: " << e.getMessage() << std::endl;
      //      std::cout << "Message: " << e.getMessage() << std::endl;
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
  std::cout << "OnCalServer::BeginTimeStamp: Setting BOR TimeStamp to " << beginTimeStamp << std::endl;
}

void OnCalServer::EndTimeStamp(const PHTimeStamp &TimeStp)
{
  endTimeStamp = TimeStp;
  std::cout << "OnCalServer::EndTimeStamp: Setting EOR TimeStamp to " << endTimeStamp << std::endl;
}

PHTimeStamp *
OnCalServer::GetLastGoodRunTS(OnCal *calibrator, const int irun)
{
  PHTimeStamp *ts = nullptr;
  if (!connectDB())
  {
    std::cout << "could not connect to " << database << std::endl;
    return ts;
  }
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  std::ostringstream subsystable;
  subsystable << "oncal" << calibrator->Name();
  cmd << "SELECT runnumber FROM " << successTable << " where runnumber < "
      << irun << " and "
      << calibrator->Name() << " > 0 order by runnumber desc limit 1";
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << subsystable.str() << " does not exist" << std::endl;
    return ts;
  }
  if (rs->next())
  {
    RunToTime *rt = RunToTime::instance();
    int oldrun = rs->getInt("runnumber");
    ts = rt->getBeginTime(oldrun);
    std::cout << "Getting previous good run, current run: " << irun
              << ", previous good run: " << oldrun
              << " began ";
    ts->print();
    std::cout << std::endl;
  }
  else
  {
    std::cout << PHWHERE << " No previous good run found for run " << irun << std::endl;
  }
  delete rs;
  return ts;
}

int OnCalServer::SyncCalibTimeStampsToOnCal(const OnCal *calibrator, const int commit)
{
  std::vector<std::string> caltab;
  calibrator->GetPdbCalTables(caltab);
  std::vector<std::string>::const_iterator iter;
  for (iter = caltab.begin(); iter != caltab.end(); ++iter)
  {
    std::cout << "dealing with table: " << *iter << std::endl;
    SyncCalibTimeStampsToOnCal(calibrator, *iter, commit);
  }
  return 0;
}

int OnCalServer::SyncCalibTimeStampsToOnCal(const OnCal *calibrator, const std::string &table, const int commit)
{
  std::string name = calibrator->Name();
  odbc::Connection *con = nullptr;
  odbc::Connection *concalib = nullptr;
  std::ostringstream cmd;
  try
  {
    con = odbc::DriverManager::getConnection(database.c_str(), "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot connect to " << database.c_str() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }
  try
  {
    concalib = odbc::DriverManager::getConnection("oncal", "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot connect to "
              << "oncal" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }
  odbc::Statement *stmt = nullptr;
  try
  {
    stmt = con->createStatement();
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  odbc::PreparedStatement *stmt1 = nullptr;
  odbc::ResultSet *rs1 = nullptr;
  try
  {
    cmd.str("");
    cmd << "SELECT * from " << table << " where startvaltime = ?";
    stmt1 = concalib->prepareStatement(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  odbc::PreparedStatement *stmtupd = nullptr;
  try
  {
    cmd.str("");
    cmd << "update " << table << " set endvaltime = ? where startvaltime = ?";
    stmtupd = concalib->prepareStatement(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  cmd.str("");
  cmd << "select * from "
      << successTable
      << " where "
      << name
      << " > 0";
  //      << " > 0 and runnumber < 150000";
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
    return -1;
  }
  while (rs->next())
  {
    int run = rs->getInt("runnumber");
    int startticks = rs->getLong("startvaltime");
    int endticks = rs->getLong("endvaltime");
    // int status = rs->getInt(name);
    //       std::cout << "run: " << run
    // 	   << ", status: " << status
    // 	   << ", startticks: " << startticks
    // 	   << ", endticks: " << endticks << std::endl;
    stmt1->setInt(1, startticks);
    try
    {
      rs1 = stmt1->executeQuery();
    }
    catch (odbc::SQLException &e)
    {
      std::cout << "Message: " << e.getMessage() << std::endl;
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
          std::cout << "endvaltime problem with run " << run << std::endl;
          std::cout << "endvaltime from oncal_status: " << endticks << std::endl;
          std::cout << "startvaltime from oncal_status: " << startticks << std::endl;
          std::cout << "endvaltime from calibrations DB: " << rs1->getInt("endvaltime") << std::endl;
          if (endticks < rs1->getInt("endvaltime"))
          {
            std::cout << "ENDTICKS smaller CALIB" << std::endl;
            //                      return -1;
          }
        }
        isproblem = 1;
      }
      else
      {
        if (isproblem)
        {
          std::cout << "endvaltime changes, check run " << run << std::endl;
          //                  return -1;
        }
      }
      // 	   std::cout << "starttime: " << rs1->getInt("startvaltime") << std::endl;
      // 	   std::cout << "endtime: " << rs1->getInt("endvaltime") << std::endl;
      ionce++;
    }
    if (isproblem)
    {
      std::cout << "Adjusting run " << run << std::endl;
      std::cout << "changing endvaltime from " << calibendval
                << " to " << endticks << std::endl;
      if (commit)
      {
        stmtupd->setInt(1, endticks);
        stmtupd->setInt(2, startticks);
        stmtupd->executeUpdate();
      }
    }
    if (!ionce)
    {
      std::cout << "Run " << run << " not found" << std::endl;
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
  std::ostringstream cmd;
  try
  {
    con = odbc::DriverManager::getConnection(database.c_str(), "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot connect to " << database.c_str() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }
  odbc::Statement *stmt = nullptr;
  try
  {
    stmt = con->createStatement();
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  odbc::PreparedStatement *stmtupd = nullptr;
  try
  {
    cmd.str("");
    cmd << "UPDATE oncal_status set endvaltime = ? where runnumber = ?";
    stmtupd = con->prepareStatement(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  cmd.str("");
  cmd << "select * from "
      << successTable;  //<< " where runnumber > 160000";
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Message: " << e.getMessage() << std::endl;
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
      std::cout << "Run " << run
                << ": Start mismatch, oncal: " << startticks
                << ", rt: " << rtstartticks << std::endl;
    }
    if (rtendticks != endticks)
    {
      // exclude starttime=endtime in runtotime (some crashed calibrations can do this)
      // in this case the calibration adds 1 sec to starttime
      if (rtstartticks != rtendticks)
      {
        std::cout << "Run " << run
                  << ": End mismatch, oncal: " << endticks
                  << ", rt: " << rtendticks << std::endl;
        if (endticks > rtendticks)
        {
          std::cout << "BAD: endticks: " << endticks
                    << ", rtendticks: " << rtendticks
                    << std::endl;
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
          std::cout << "Run " << run
                    << ": Start/End mismatch, Start: " << startticks
                    << ", End: " << endticks << std::endl;
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
            std::cout << "run " << run << " was twiddled by OnCal" << std::endl;
          }
        }
      }
    }
    //       std::cout << "run: " << run
    // 	   << ", status: " << status
    // 	   << ", startticks: " << startticks
    // 	   << ", endticks: " << endticks << std::endl;
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

int OnCalServer::CreateCalibration(OnCal *calibrator, const int myrunnumber, const std::string &what, const int commit)
{
  int iret = -1;
  runNum = myrunnumber;
  SetBorTime(myrunnumber);
  SetEorTime(myrunnumber);
  if (!connectDB())
  {
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  add_calibrator_to_statustable(calibrator->Name());
  std::string table = "OnCal";
  table += calibrator->Name();
  check_create_subsystable(table);
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;

  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << myrunnumber;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
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
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE << " Exception caught, Message: "
              << e.getMessage() << std::endl;
    return -1;
  }
  if (rs->next() || testmode)
  {
    PdbBankManager *bankManager = PdbBankManager::instance();
    PdbApplication *application = bankManager->getApplication();
    application->setDBName("oncal");
    std::string tablecomment = "Subsytem provided";
    iret = calibrator->CreateCalibration(runnumber, what, tablecomment, commit);
    if (!iret)
    {
      std::cout << "Comment: " << tablecomment << std::endl;
      std::cout << "updating oncal status tables for " << runnumber << std::endl;
      if (commit)
      {
        CreateCalibrationUpdateStatus(calibrator, table, tablecomment, OnCalDBCodes::SUBSYSTEM);
      }
    }
    else
    {
      std::cout << "Calibratior " << calibrator->Name() << " for run " << runnumber << " failed" << std::endl;
      if (commit)
      {
        CreateCalibrationUpdateStatus(calibrator, table, tablecomment, OnCalDBCodes::FAILED);
      }
    }
  }
  else
  {
    std::cout << PHWHERE << " Run " << runnumber << " is already successfully calibrated for "
              << calibrator->Name() << std::endl;
  }
  return iret;
}

void OnCalServer::CreateCalibrationUpdateStatus(OnCal *calibrator, const std::string &table, const std::string &tablecomment, const int dbcode)
{
  updateDB(successTable.c_str(), calibrator->Name(), dbcode);
  insertRunNumInDB(table, RunNumber());
  updateDB(table, "comment", tablecomment, RunNumber(), true);
  std::ostringstream stringarg;
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
  std::vector<std::string> flist = calibrator->GetLocalFileList();
  if (flist.size())
  {
    std::string filelist = "";
    for (const std::string &infile : flist)
    {
      filelist += infile;
      filelist += " ";
    }
    filelist.pop_back();  // strip empty space at end from loop
    std::cout << "FileList: " << filelist << std::endl;
    updateDB(table, "files", filelist, RunNumber());
  }
  return;
}

int OnCalServer::CopySnglTable(const std::string &pdbclass, const std::string &tablename, const int bankid, const int FromRun, const int ToRun, const int commit)
{
  int iret = CopySnglTableNewBankId(pdbclass, tablename, bankid, bankid, FromRun, ToRun, commit);
  return iret;
}

int OnCalServer::CopySnglTableNewBankId(const std::string &pdbclass, const std::string &tablename, const int bankid, const int Tobankid, const int FromRun, const int ToRun, const int commit)
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
  std::ostringstream newdesc;
  newdesc << "copied from run " << FromRun;
  if (pdbBanknew)
  {
    pdbBanknew->setLength(pdbBank->getLength());
    if (Verbosity() > 0)
    {
      for (unsigned int i = 0; i < pdbBank->getLength(); i++)
      {
        std::cout << "orig: " << std::endl;
        pdbBank->printEntry(i);
        std::cout << "new: " << std::endl;
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
      std::cout << "StartValTime: orig " << pdbBank->getStartValTime()
                << ", new: " << pdbBanknew->getStartValTime() << std::endl;
      std::cout << "EndValTime: orig " << pdbBank->getEndValTime()
                << ", new: " << pdbBanknew->getEndValTime() << std::endl;
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
    std::cout << PHWHERE << "Unknown Run " << irun << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;

  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber < "
      << irun << " and "
      << calibrator->Name() << " = 1 order by runnumber desc limit 1";
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
    return -1;
  }
  int prevrun = -1;
  unsigned int prevend = 0;
  if (rs->next())
  {
    prevrun = rs->getInt("runnumber");
    unsigned int prevstart = rs->getLong("startvaltime");
    prevend = rs->getLong("endvaltime");
    std::cout << "previous run: " << prevrun
              << ", start: " << prevstart
              << ", end: " << prevend
              << std::endl;
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << " No previous good run found for run " << irun << std::endl;
    }
  }
  delete rs;
  closestrun = prevrun;
  if (previous == fetchrun::PREVIOUS)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Closest previous run is " << closestrun << std::endl;
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
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
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
      std::cout << "next run: " << nextrun
                << ", start: " << nextstart
                << ", end: " << nextend
                << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << " No next good run found for run " << irun << std::endl;
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
    std::cout << "diff prev: " << tdiffprev
              << ", next: " << tdiffnext
              << std::endl;
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
    std::cout << "closest run: " << closestrun << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  runNum = runno;
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;

  cmd << "SELECT runnumber FROM "
      << successTable << " where runnumber = "
      << runno;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
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
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE << " Exception caught, Message: "
              << e.getMessage() << std::endl;
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
        std::cout << "ClosestGoodRun returned bad runnumber: " << FromRun << std::endl;
        return -1;
      }
    }
    std::cout << "Going to copy calibration for run " << runno
              << " from run " << FromRun << std::endl;

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
      std::string table = "OnCal";
      table += calibrator->Name();
      std::ostringstream comment;
      comment << " CopiedRun(" << FromRun << ")";
      std::cout << "updating oncal status tables for " << runno << std::endl;
      if (commit)
      {
        updateDB(successTable.c_str(), calibrator->Name(), newstatus);
        insertRunNumInDB(table, runNum);
        updateDB(table, "comment", comment.str(), runNum, true);
        updateDB(table, "committed", true);
      }
    }
  }
  else
  {
    std::cout << "Run " << runno
              << " has a good calibrations, doing nothing" << std::endl;
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
    std::cout << PHWHERE << "Cannot get begin time for run " << runno << std::endl;
    std::cout << "Exiting" << std::endl;
    exit(1);
  }
  BeginTimeStamp(*BorTimeStp);

  // enter begin run timestamp into rc flags
  PHTimeStamp BeginRunTimeStamp(*BorTimeStp);
  //  rc->set_TimeStamp(BeginRunTimeStamp);
  std::cout << "OnCalServer::SetBorTime from RunToTime was found for run : " << runno << " to ";
  BeginRunTimeStamp.print();
  std::cout << std::endl;

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
  std::cout << "OnCalServer::SetEorTime: setting eor time to ";
  EorTimeStp->print();
  std::cout << std::endl;
  delete EorTimeStp;
  return 0;
}

int OnCalServer::GetRunTimeTicks(const int runno, time_t &borticks, time_t &eorticks)
{
  RunToTime *runTime = RunToTime::instance();
  PHTimeStamp *TimeStp(runTime->getBeginTime(runno));
  if (!TimeStp)
  {
    std::cout << PHWHERE << "Cannot get begin time for run " << runno << std::endl;
    std::cout << "Exiting" << std::endl;
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
  std::map<std::string, std::set<SubsysReco *> >::iterator iter;
  if (check_calibrator_in_statustable(calibratorname))
  {
    std::cout << PHWHERE << " the calibrator " << calibratorname << " is unknown to me" << std::endl;
    return -1;
  }
  iter = requiredCalibrators.find(calibratorname);
  if (iter != requiredCalibrators.end())
  {
    iter->second.insert(reco);
  }
  else
  {
    std::set<SubsysReco *> subsys;
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
    std::cout << PHWHERE << "Unknown Run " << irun << std::endl;
    return -1;
  }
  if (requiredCalibrators.size() == 0)
  {
    std::cout << PHWHERE << "No required calibrations given" << std::endl;
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
    std::cout << "could not connect to " << database << std::endl;
    return -1;
  }
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;
  std::map<std::string, std::set<SubsysReco *> >::const_iterator iter;
  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT runnumber,startvaltime,endvaltime FROM "
      << successTable << " where runnumber <= "
      << irun;
  for (iter = requiredCalibrators.begin(); iter != requiredCalibrators.end(); ++iter)
  {
    cmd << " and " << iter->first << " > 0 ";
  }

  cmd << " order by runnumber desc limit 1";
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
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
      std::cout << "previous run: " << prevrun
                << ", start: " << prevstart
                << ", end: " << prevend
                << std::endl;
    }
  }
  else
  {
    std::cout << PHWHERE << " No previous good run found for run " << irun << std::endl;
  }
  delete rs;
  // if the current run fullfills requirements return immediately
  if (prevrun == irun)
  {
    std::cout << "closest run with required calibs is current run: " << irun << std::endl;
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
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
    return -1;
  }
  int nextrun = 0;
  unsigned int nextstart = 0;
  if (rs->next())
  {
    nextrun = rs->getInt("runnumber");
    nextstart = rs->getLong("startvaltime");
    unsigned int nextend = rs->getLong("endvaltime");
    std::cout << "next run: " << nextrun
              << ", start: " << nextstart
              << ", end: " << nextend
              << std::endl;
  }
  else
  {
    std::cout << PHWHERE << " No next good run found for run " << irun << std::endl;
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
  std::cout << "closest run with required calibs: " << closestrun << std::endl;
  return closestrun;
}

int OnCalServer::FillRunListFromFileList()
{
  for (Fun4AllSyncManager *sync : SyncManagers)
  {
    for (Fun4AllInputManager *inmgr : sync->GetInputManagers())
    {
      for (const std::string &infile : inmgr->GetFileList())
      {
        std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(infile);
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
  std::string table = "OnCalRichCal";
  check_create_subsystable(table);
  GetRunTimeTicks(firstrun, beginticks, dummy);
  GetRunTimeTicks(lastrun, dummy, endticks);
  std::ostringstream stringarg;
  stringarg << OnCalDBCodes::COVERED;
  //  std::set<int>::const_iterator runiter;
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
    std::string tablename = "calibrichadc";
    odbc::Connection *con = 0;
    std::ostringstream cmd;
    try
    {
    con = odbc::DriverManager::getConnection("oncal", "phnxrc", "");
    }
    catch (odbc::SQLException& e)
    {
    std::cout << "Cannot connect to " << database.c_str() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
    }
    odbc::Statement *stmt = 0;
    odbc::Statement *stmtup = 0;
    try
    {
    stmt = con->createStatement();
    stmtup = con->createStatement();
    }
    catch (odbc::SQLException& e)
    {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
    }

    odbc::ResultSet *rs1 = 0;
    cmd.str("");
    cmd << "SELECT  endvaltime from " << tablename
    << " where bankid = 1 and startvaltime = " << beginticks;
    std::cout << "sql cmd: " << cmd.str() << std::endl;
    try
    {
    rs1 = stmt->executeQuery(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
    }
    if (rs1->next())
    {
    std::cout << "Endcaltime: " << rs1->getInt("endvaltime") << std::endl;
    std::cout << "future endvaltime: " << endticks << std::endl;
    cmd.str("");
    cmd << "Update " << tablename
    << " set endvaltime = " << endticks
    << " where bankid = 1 and startvaltime = "
    << beginticks;
    stmtup->executeUpdate(cmd.str());

    }
    else
    {
    std::cout << "Could not find startvaltime " << beginticks
    << "from run " << firstrun << std::endl;
    }

  */

  return 0;
}

int OnCalServer::GetCalibStatus(const std::string &calibname, const int runno)
{
  int iret = -3;
  if (!connectDB())
  {
    std::cout << "could not connect to " << database << std::endl;
    return -4;
  }
  odbc::Statement *stmt = DBconnection->createStatement();
  std::ostringstream cmd;

  // look only for runs which were actually successfully calibrated (status = 1)
  cmd << "SELECT " << calibname << " FROM "
      << successTable << " where runnumber = "
      << runno;
  std::cout << "exec " << cmd.str() << std::endl;
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Table " << successTable << " does not exist" << std::endl;
    return -5;
  }
  if (rs->next())
  {
    iret = rs->getInt(calibname.c_str());
  }
  else
  {
    std::cout << PHWHERE << " No calib status for " << calibname
              << " for " << runno << std::endl;
  }
  delete rs;
  return iret;
}

void OnCalServer::TestMode(const int i)
{
  const char *logname = getenv("LOGNAME");
  if (logname)
  {
    if (strcmp(logname, "sphnxpro") == 0 || strcmp(logname, "anatrain") == 0)
    {
      std::cout << "phnxcal,anatrain account is not allowed to run in testmode" << std::endl;
    }
    else
    {
      testmode = i;
    }
  }
  else
  {
    std::cout << "could not get account via env var LOGNAME, not setting testmode" << std::endl;
  }
  return;
}

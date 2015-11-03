// $Id: PgPostBankBackupLog.cc,v 1.2 2014/05/19 17:06:23 jinhuang Exp $

/*!
 * \file PgPostBankBackupLog.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2014/05/19 17:06:23 $
 */

#include "PgPostBankManager.hh"
#include "PgPostHelper.hh"
#include "PgPostCalBankIterator.hh"
#include "PgPostBankWrapper.hh"
#include "PgPostBankWrapper2.hh"
#include "PgPostApplication.hh"
#include "PgPostCalBank.hh"

#include <PdbBankID.hh>
#include <PdbBankID2.hh>
#include <PHString.h>
#include <PHPointerList.h>
#include <RDBC/TSQL.h>
#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLResultSetMetaData.h>
#include <RDBC/TSQLPreparedStatement.h>
#include <RDBC/TSQLDatabaseMetaData.h>
#include <PdbBankList.hh>
#include <PdbCalBank.hh>
#include <PdbClassMap.hh>
#include <PHString.h>
#include <PHPointerList.h>
#include <RunToTimePg.hh>
#include <PHTimeServer.h>
#include <PdbBankManagerFactory.hh>

#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TBufferFile.h>

#include <ctime>
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <algorithm>

using namespace std;

#include "PgPostBankBackupLog.hh"

PgPostBankBackupLog::TSQLConnection_PTR PgPostBankBackupLog::con(
    static_cast<TSQLConnection *>(NULL));

PgPostBankBackupLog::PgPostBankBackupLog(const std::string & TableName,
    const std::string & Tag) :
    verbosity(1), tablename(TableName), tag(Tag), pstmt(NULL)
{

}

PgPostBankBackupLog::~PgPostBankBackupLog()
{

}

void
PgPostBankBackupLog::Init()
{

  try
    {
      if (!con)
        {
          const TString s_con = "dsn=calibrations_backup_log; uid=phnxrc; pwd= ";

          if (Verbosity() >= 1)
            cout << "PgPostBankBackupLog::Init - connect to " << s_con << endl;

          con = gSQLDriverManager->GetConnection(s_con);
          if (!con)
            {
              cout << "PgPostBankBackupLog::Init - Error - "
                  << "cannot init connection to " << s_con << endl;
              exit(1);
            }
        }

      if (!pstmt)
        {
          assert(con);

          ostringstream sqlcmd;
          sqlcmd << "insert into calib_log (tablename,rid,ops,tag) values ('"
              << tablename << "',?,?,'" << tag << "');";

          if (Verbosity() >= 1)
            cout << "PgPostBankBackupLog::Init - make TSQLPreparedStatement "
                << sqlcmd.str() << endl;
          pstmt = con->PrepareStatement(sqlcmd.str().c_str());

          if (!pstmt)
            {
              cout << "PgPostBankBackupLog::Init - Error - "
                  << "cannot prepare statement " << sqlcmd.str() << endl;
              exit(1);
            }
        }
    }
  catch (std::exception& e)
    {
      cout << "PgPostBankBackupLog::Init - Error - " << "Initialization error "
          << e.what() << endl;
      exit(1);
    }
}

void
PgPostBankBackupLog::Log(const int rid, const enu_ops ops)
{
  Init();

  if (Verbosity() >= 2)
    {
      cout << "PgPostBankBackupLog::Log - log rid = " << rid << ", ops = "
          << ops << " for table " << tablename << " and tag " << tag << endl;
    }

  assert(con);
  assert(pstmt);

  pstmt->SetInt(1, rid);
  pstmt->SetInt(2, (int) ops);
  try
    {
      const int res = pstmt->ExecuteUpdate();

      if (res == 0)
        {
          cout << "PgPostBankBackupLog::Log - Error - "
              << "DATABASE: commit failed. "
              << "Make sure you commit to the writable database " << endl;
          exit(1);
        }

      con->Commit();
    }
  catch (TSQLException& e)
    {
      cout << "PgPostBankBackupLog::Log - Error - "
          << " Exception caught during connection->commit()" << endl;
      cout << e.GetMessage() << endl;
      exit(1);
    }

}

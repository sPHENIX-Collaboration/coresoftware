#include "PgPostBankBackupLog.h"
#include "PgPostApplication.h"
#include "PgPostBankManager.h"
#include "PgPostBankWrapper.h"
#include "PgPostCalBank.h"
#include "PgPostCalBankIterator.h"
#include "RunToTimePg.h"

#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankList.h>
#include <pdbcalbase/PdbBankManagerFactory.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbClassMap.h>

#include <phool/PHPointerList.h>
#include <phool/PHTimeServer.h>

#include <RDBC/TSQL.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLDatabaseMetaData.h>
#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLPreparedStatement.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLResultSetMetaData.h>

#include <TBufferFile.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

using namespace std;

PgPostBankBackupLog::TSQLConnection_PTR PgPostBankBackupLog::con(
    static_cast<TSQLConnection*>(nullptr));

PgPostBankBackupLog::PgPostBankBackupLog(const std::string& TableName,
                                         const std::string& Tag)
  : verbosity(1)
  , tablename(TableName)
  , tag(Tag)
  , pstmt(nullptr)
{
}

PgPostBankBackupLog::~PgPostBankBackupLog()
{
}

void PgPostBankBackupLog::Init()
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
    cout << "PgPostBankBackupLog::Init - Error - "
         << "Initialization error "
         << e.what() << endl;
    exit(1);
  }
}

void PgPostBankBackupLog::Log(const int rid, const enu_ops ops)
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

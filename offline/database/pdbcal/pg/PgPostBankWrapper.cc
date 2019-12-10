#include "PgPostBankWrapper.h"
#include "PgPostApplication.h"
#include "PgPostBankWrapperManager.h"

#include <pdbcalbase/PdbBankID.h>        // for PdbBankID
#include <pdbcalbase/PdbCalBank.h>       // for PdbCalBank

#include <phool/phool.h>
#include <phool/PHTimeStamp.h>

#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLPreparedStatement.h>

#include <unistd.h>  // for sleep
#include <cstdlib>
#include <iostream>
#include <sstream>


using namespace std;

PgPostBankWrapper::PgPostBankWrapper()
  : bank(nullptr)
{
}

PgPostBankWrapper::PgPostBankWrapper(PdbCalBank *b)
{
  bank = b;
  PgPostBankWrapperManager::instance().registerWrapper(this);
}

PgPostBankWrapper::~PgPostBankWrapper()
{
  delete bank;
  PgPostBankWrapperManager::instance().unregisterWrapper(this);
}

void PgPostBankWrapper::printHeader() const
{
  bankID.print();
  cout << "Insert   : " << insertTime << endl;
  cout << "StartVal : " << startValTime << endl;
  cout << "EndVal   : " << endValTime << endl;
  cout << "Calibration Description : " << description << endl;
  cout << "User Name = " << userName << endl;
}

bool PgPostBankWrapper::commit()
{
  if (bank)
  {
    PgPostApplication *ap = PgPostApplication::instance();
    TSQLConnection *con = ap->getConnection();

    ostringstream sqlcmd;
    sqlcmd << "insert into " << ((*this).getTableName())
           << " values(?,?,?,?,?,?,?);";
    cout << "query: " << sqlcmd.str() << endl;
    TSQLPreparedStatement *pstmt = con->PrepareStatement(sqlcmd.str().c_str());
    pstmt->SetInt(1, ((*this).getBankID()).getInternalValue());
    pstmt->SetLong(2, ((*this).getInsertTime()).getTics());
    pstmt->SetLong(3, ((*this).getStartValTime()).getTics());
    pstmt->SetLong(4, ((*this).getEndValTime()).getTics());
    pstmt->SetString(5, ((*this).getDescription()));
    pstmt->SetString(6, ((*this).getUserName()));
    pstmt->SetObject(7, this);
    int res = 0;
    res = pstmt->ExecuteUpdate();
    con->Commit();
    if (res == 0)
    {
      cout << PHWHERE << "DATABASE: commit to " << tableName << " failed" << endl;
      cout << "Make sure you commit to the master database on phnxdb2.phenix.bnl.gov" << endl;
      cout << "and that " << tableName << " exists in the master database" << endl;
      exit(1);
    }
    cout << "Committed " << res << " row, sleeping 1 sec to make sure we have distinct insert times" << endl;
    // entries are distinguished by insert time. If one commits a lot at once one has multiple
    // entries for the same insert time (we can commit ~4/sec). This sleep 1 second just forces
    // the insert time to be unique.
    sleep(1);
    return res;
  }
  else
  {
    cout << "Bank is a nullptr pointer" << endl;
    return 0;
  }
}

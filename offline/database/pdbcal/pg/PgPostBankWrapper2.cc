#include "PgPostBankWrapper2.h"

#include <phool/phool.h>

#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLResultSetMetaData.h>
#include <RDBC/TSQLPreparedStatement.h>
#include <PgPostApplication.h>
#include <PgPostBankWrapperManager.h>

#include <RDBC/TSQLDatabaseMetaData.h>

#include <cstdlib>
#include <iostream>
#include <unistd.h> // for sleep
#include <sstream>

using namespace std;

PgPostBankWrapper2::PgPostBankWrapper2()
{
  bank = 0;
}

PgPostBankWrapper2::PgPostBankWrapper2(PdbCalBank *b)
{
  bank = b;
  PgPostBankWrapperManager::instance().registerWrapper2(this);
}

PgPostBankWrapper2::~PgPostBankWrapper2()
{
  delete bank;
  PgPostBankWrapperManager::instance().unregisterWrapper2(this);
}

void
PgPostBankWrapper2::printHeader() const
{
  bankID2.print();
  cout << "Insert   : " << insertTime << endl;
  cout << "StartVal : " << startValTime << endl;
  cout << "EndVal   : " << endValTime << endl;
  cout << "Calibration Description : " << description << endl;
  cout << "User Name = " << userName << endl;
}

bool PgPostBankWrapper2::commit()
{
  if (bank)
    {
      PgPostApplication *ap = PgPostApplication::instance();
      TSQLConnection *con = ap->getConnection();
      
      ostringstream sqlcmd;
      sqlcmd << "insert into " << ((*this).getTableName())
	<< " values(?,?,?,?,?,?,?);";
      cout << "query: " << sqlcmd.str() << endl;
      TSQLPreparedStatement* pstmt = con->PrepareStatement(sqlcmd.str().c_str());
      pstmt->SetInt(1, ((*this).getBankID2()).getInternalValue());
      pstmt->SetLong(2, ((*this).getInsertTime()).getTics());
      pstmt->SetLong(3, ((*this).getStartValTime()).getTics());
      pstmt->SetLong(4, ((*this).getEndValTime()).getTics());
      pstmt->SetString(5, ((*this).getDescription()));
      pstmt->SetString(6, ((*this).getUserName()));
      pstmt->SetObject(7, this);
      int res = 0;
      res = pstmt->ExecuteUpdate();
    try
      {
  con->Commit();
      }
    catch (TSQLException& e)
    {
      cout << PHWHERE
	   << " Exception caught during connection->commit()" << endl;
      cout << e.GetMessage() << endl;
      return 1;
    }
      if (res == 0){
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
      cout << "Bank is a NULL pointer" << endl;
      return 0;
    }
}

bool PgPostBankWrapper2::commit_rid(int rid,long it, long st,long et)
{
  if (bank)
    {
      PgPostApplication *ap = PgPostApplication::instance();
      TSQLConnection *con = ap->getConnection();
      
      ostringstream sqlcmd;
      sqlcmd << "insert into " << ((*this).getTableName())
	<< " values(?,?,?,?,?,?,?,?);";
      cout << "query: " << sqlcmd.str() << endl;
      TSQLPreparedStatement* pstmt = con->PrepareStatement(sqlcmd.str().c_str());
      pstmt->SetInt(1, ((*this).getBankID2()).getInternalValue());
      pstmt->SetLong(2, it);
      pstmt->SetLong(3, st);
      pstmt->SetLong(4, et);
      pstmt->SetString(5, ((*this).getDescription()));
      pstmt->SetString(6, ((*this).getUserName()));
      pstmt->SetObject(7, this);
      pstmt->SetInt(8,rid);
      int res = 0;
      try{
      res = pstmt->ExecuteUpdate();
      }
catch (TSQLException& e)
	{
	  cout << PHWHERE
	       << " Exception caught during connection->commit()" << endl;
	  cout << e.GetMessage() << endl;
	  return 1;
	}
      try{
	con->Commit();
      }
      catch (TSQLException& e)
	{
	  cout << PHWHERE
	       << " Exception caught during connection->commit()" << endl;
	  cout << e.GetMessage() << endl;
	  return 1;
	}
      if (res == 0){
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
      cout << "Bank is a NULL pointer" << endl;
      return 0;
    }
}

ClassImp(PgPostBankWrapper2)


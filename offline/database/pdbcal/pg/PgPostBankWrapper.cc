#include "PgPostBankWrapper.hh"
#include "PgPostApplication.hh"
#include "PgPostBankWrapperManager.hh"

#include <phool.h>

#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLResultSetMetaData.h>
#include <RDBC/TSQLPreparedStatement.h>

#include <RDBC/TSQLDatabaseMetaData.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h> // for sleep

ClassImp(PgPostBankWrapper)

using namespace std;

PgPostBankWrapper::PgPostBankWrapper():
  bank(NULL)
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

void
PgPostBankWrapper::printHeader() const
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
      sqlcmd << "insert into " << ((*this).getTableName()).getString()
	<< " values(?,?,?,?,?,?,?);";
      cout << "query: " << sqlcmd.str() << endl;
      TSQLPreparedStatement* pstmt = con->PrepareStatement(sqlcmd.str().c_str());
      pstmt->SetInt(1, ((*this).getBankID()).getInternalValue());
      pstmt->SetLong(2, ((*this).getInsertTime()).getTics());
      pstmt->SetLong(3, ((*this).getStartValTime()).getTics());
      pstmt->SetLong(4, ((*this).getEndValTime()).getTics());
      pstmt->SetString(5, ((*this).getDescription()).getString());
      pstmt->SetString(6, ((*this).getUserName()).getString());
      pstmt->SetObject(7, this);
      int res = 0;
      res = pstmt->ExecuteUpdate();
      con->Commit();
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

void
PgPostBankWrapper::setDescription(const PHString & val)
{
  // strncpy does not append \0 if number of characters in input string
  // exceed number of chars in output string
  // so we set the last char to \0 by hand
  strncpy(description, val.getString(),sizeof(description)-1);
  description[sizeof(description)-1] = '\0';
  if (strlen(val.getString()) > strlen(description))
    {
      cout << "description string length " << strlen(val.getString())
	   << " exceeds maximum length of " << sizeof(description)
	   << " description used in DB: " << endl << description
	   << endl;
    }
}

void
PgPostBankWrapper::setUserName(const PHString & val)
{
  // strncpy does not append \0 if number of characters in input string
  // exceed number of chars in output string
  // so we set the last char to \0 by hand
  strncpy(userName, val.getString(),sizeof(userName)-1);
  userName[sizeof(userName)-1] = '\0';
  if (strlen(val.getString()) > strlen(userName))
    {
      cout << "userName string length " << strlen(val.getString())
	   << " exceeds maximum length of " << sizeof(userName)
	   << " userName used in DB: " << endl << userName
	   << endl;
    }
}

void
PgPostBankWrapper::setTableName(const PHString & val)
{
  if (strlen(val.getString()) > sizeof(tableName))
    {
      cout << "length of tablename " << val
	   << " exceeds max length " << sizeof(tableName)
	   << " fatal, exiting now" << endl;
      exit(1);
    }
  // the above if should take care of overflows and the following use of strncpy
  // is not needed. I just want to avoid the use of strcpy
  strncpy(tableName, val.getString(),sizeof(tableName)-1);
  tableName[sizeof(tableName)-1] = '\0';
}

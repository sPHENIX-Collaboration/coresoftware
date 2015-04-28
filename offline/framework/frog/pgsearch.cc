#include "pgsearch.h"

#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/setup.h>
#include <odbc++/types.h>
#include <odbc++/errorhandler.h>
#include <sql.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/preparedstatement.h>
#include <odbc++/databasemetadata.h>

#include <cstring>
#include <iostream>
#include <string>

using namespace odbc;
using namespace std;

void
pgsearch::search(const string &lname, string &cp)
{
  Connection* con = NULL;
  Statement* stmt;
  ResultSet* rs;
  cp = lname; // if things fail, return input string
  try
    {
      con = DriverManager::getConnection("FileCatalog", "argouser", "Brass_Ring");
    }
  catch (SQLException& e)
    {
      cout << PHWHERE
	   << " Exception caught during DriverManager::getConnection" << endl;
      cout << "Message: " << e.getMessage() << endl;
      return ;
    }


  string temp = lname;
  string mys = "SELECT * from files where lfn='" + temp + "' and full_host_name <> 'hpss'";

  stmt = con->createStatement();
  rs = stmt->executeQuery(mys.c_str());

  if (rs->next())
    {
      cp = rs->getString(3);
    }
  else
    {
      delete rs;
      //return empty str if dCache instance found, hpss path if not
      mys = "SELECT * from files where lfn='" + temp + "' and full_host_name ='hpss' and full_file_path like '/home/dcphenix%'";
      try
        {
          rs = stmt->executeQuery(mys.c_str());
        }
      catch (SQLException& e)
        {
          cout << PHWHERE
	       << " Exception caught during DriverManager::getConnection" << endl;
          cout << "Message: " << e.getMessage() << endl;
        }
      if (rs->next())
        {
          cp = "";
        }
      else {
	delete rs;
	mys = "SELECT * from files where lfn='" + temp + "'  and full_file_path not like '/home/dcphenix%'";
	try
	  {
	    rs = stmt->executeQuery(mys.c_str());
	  }
	catch (SQLException& e)
	  {
	    cout << PHWHERE
		 << " Exception caught during DriverManager::getConnection" << endl;
	    cout << "Message: " << e.getMessage() << endl;
	  }
	if (rs->next())
	  {
	    cp = rs->getString(3);
	  }
	else {
	  cp= "";
	}
      }
    }
  delete rs;
  delete con;
}

#include "dCachesearch.h"

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

#include <sys/stat.h>

#include <cstring> 
#include <iostream>
#include <string>

using namespace odbc;
using namespace std;


void
dCachesearch::search(const string &lname, string &cp)
{
  struct stat64 stbuf;
  Connection* con = 0;
  Statement* stmt;
  ResultSet* rs;
  const char* dcp;
  string dc;
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
  string mys = "SELECT * from files where lfn='" + temp + "' and full_host_name = 'hpss' and full_file_path like '/home/dcphenix/phnxreco/%'";

  stmt = con->createStatement();
  rs = stmt->executeQuery(mys.c_str());

  if (rs->next())
    {
      cp = rs->getString(3);
      cp = cp.substr(14, cp.size());
      dc = "/direct/phenix+pnfs" + cp;
      dcp = dc.c_str();
      if (stat64(dcp, &stbuf) != -1)
        {
          if ((stbuf.st_mode & S_IFMT) == S_IFREG)
            {
              //cp = "dcap://dcphenix01.rcf.bnl.gov:22125//pnfs/rcf.bnl.gov/phenix"+cp;
              cp = "dcache:/direct/phenix+pnfs" + cp;
            }
        }
      else
        {
          cp = "/home" + cp;
        }
    }
  else
    {
      cp = "";
    }
  delete rs;
  delete con;
}

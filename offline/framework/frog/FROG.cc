#include "FROG.h"

#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException

#include <boost/tokenizer.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

const char *
FROG::location(const std::string &logical_name)
{
  pfn = logical_name;
  if (logical_name.empty() || logical_name.find('/') != std::string::npos)
  {
    if (Verbosity() > 0)
    {
      if (logical_name.empty())
      {
        std::cout << "FROG: empty string as filename" << std::endl;
      }
      else if (logical_name.find('/') != std::string::npos)
      {
        std::cout << "FROG: found / in filename, assuming it contains a full path" << std::endl;
      }
    }
    return pfn.c_str();
  }
  try
  {
    char *gsearchpath_env = getenv("GSEARCHPATH");
    if (gsearchpath_env == nullptr)
    {
      return pfn.c_str();
    }
    std::string gsearchpath(gsearchpath_env);
    if (Verbosity() > 0)
    {
      std::cout << "FROG: GSEARCHPATH: " << gsearchpath << std::endl;
    }
    boost::char_separator<char> sep(":");
    boost::tokenizer<boost::char_separator<char> > tok(gsearchpath, sep);
    for (auto &iter : tok)
    {
      if (iter == "PG")
      {
        if (Verbosity() > 1)
        {
          std::cout << "Searching FileCatalog for disk resident file "
               << logical_name << std::endl;
        }
        if (PGSearch(logical_name))
        {
          if (Verbosity() > 1)
          {
            std::cout << "Found " << logical_name << " in FileCatalog, returning "
                 << pfn << std::endl;
          }
          break;
        }
      }
      else if (iter == "DCACHE")
      {
        if (Verbosity() > 1)
        {
          std::cout << "Searching FileCatalog for dCache file "
               << logical_name << std::endl;
        }
        if (dCacheSearch(logical_name))
        {
          if (Verbosity() > 1)
          {
            std::cout << "Found " << logical_name << " in dCache, returning "
                 << pfn << std::endl;
          }
          break;
        }
      }
      else if (iter == "XROOTD")
      {
        if (Verbosity() > 1)
        {
          std::cout << "Searching FileCatalog for XRootD file "
               << logical_name << std::endl;
        }
        if (XRootDSearch(logical_name))
        {
          if (Verbosity() > 1)
          {
            std::cout << "Found " << logical_name << " in XRootD, returning "
                 << pfn << std::endl;
          }
          break;
        }
      }
      else if (iter == "LUSTRE")
      {
        if (Verbosity() > 1)
        {
          std::cout << "Searching FileCatalog for Lustre file "
               << logical_name << std::endl;
        }
        if (LustreSearch(logical_name))
        {
          if (Verbosity() > 1)
          {
            std::cout << "Found " << logical_name << " in Lustre, returning "
                 << pfn << std::endl;
          }
          break;
        }
      }
      else if (iter == "MINIO")
      {
        if (Verbosity() > 1)
        {
          std::cout << "Searching FileCatalog for Lustre file via MinIO "
               << logical_name << std::endl;
        }
        if (MinIOSearch(logical_name))
        {
          if (Verbosity() > 1)
          {
            std::cout << "Found " << logical_name << " in Lustre, returning MinIO URL "
                 << pfn << std::endl;
          }
          break;
        }
      }
      else  // assuming this is a file path
      {
        if (Verbosity() > 0)
        {
          std::cout << "Trying path " << iter << std::endl;
        }
        std::string fullfile = iter + "/" + logical_name;
        if (localSearch(fullfile))
        {
          break;
        }
      }
    }
  }
  catch (...)
  {
    if (Verbosity() > 0)
    {
      std::cout << "FROG: GSEARCHPATH not set " << std::endl;
    }
  }
  Disconnect();
  return pfn.c_str();
}

bool FROG::localSearch(const std::string &logical_name)
{
  if (std::ifstream(logical_name))
  {
    pfn = logical_name;
    return true;
  }
  return false;
}

bool FROG::GetConnection()
{
  if (m_OdbcConnection)
  {
    return true;
  }
  int icount = 0;
  do
  {
    try
    {
      m_OdbcConnection = odbc::DriverManager::getConnection("FileCatalog", "argouser", "Brass_Ring");
      return true;
    }
    catch (odbc::SQLException &e)
    {
      std::cout << PHWHERE
           << " Exception caught during DriverManager::getConnection" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
    icount++;
    std::this_thread::sleep_for(std::chrono::seconds(30));  // sleep 30 seconds before retry
  } while (icount < 5);
  return false;
}

void FROG::Disconnect()
{
  delete m_OdbcConnection;
  m_OdbcConnection = nullptr;
}

bool FROG::PGSearch(const std::string &lname)
{
  bool bret = false;
  if (!GetConnection())
  {
    return bret;
  }
  std::string sqlquery = "SELECT full_file_path from files where lfn='" + lname + "' and full_host_name <> 'hpss' and full_host_name <> 'dcache' and full_host_name <> 'lustre'";

  if (Verbosity() > 1)
  {
    std::cout << "sql query:" << std::endl
	      << sqlquery << std::endl;
  }
  odbc::Statement *stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet *rs = stmt->executeQuery(sqlquery);

  if (rs->next())
  {
    pfn = rs->getString(1);
    bret = true;
  }
  delete rs;
  delete stmt;
  return bret;
}

bool FROG::dCacheSearch(const std::string &lname)
{
  bool bret = false;
  if (!GetConnection())
  {
    return bret;
  }
  std::string sqlquery = "SELECT full_file_path from files where lfn='" + lname + "' and full_host_name = 'dcache'";

  if (Verbosity() > 1)
  {
    std::cout << "sql query:" << std::endl
	      << sqlquery << std::endl;
  }
  odbc::Statement *stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet *rs = stmt->executeQuery(sqlquery);

  if (rs->next())
  {
    std::string dcachefile = rs->getString(1);
    if (std::ifstream(dcachefile))
    {
      pfn = "dcache:" + dcachefile;
      bret = true;
    }
  }
  delete rs;
  delete stmt;
  return bret;
}

bool FROG::XRootDSearch(const std::string &lname)
{
  bool bret = false;
  if (!GetConnection())
  {
    return bret;
  }
  std::string sqlquery = "SELECT full_file_path from files where lfn='" + lname + "' and full_host_name = 'lustre'";
  if (Verbosity() > 1)
  {
    std::cout << "sql query:" << std::endl
	      << sqlquery << std::endl;
  }
  odbc::Statement *stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet *rs = stmt->executeQuery(sqlquery);

  if (rs->next())
  {
    std::string xrootdfile = rs->getString(1);
    pfn = "root://xrdsphenix.rcf.bnl.gov/"  + xrootdfile;
    bret = true;
  }
  delete rs;
  delete stmt;
  return bret;
}

bool FROG::LustreSearch(const std::string &lname)
{
  bool bret = false;
  if (!GetConnection())
  {
    return bret;
  }
  std::string sqlquery = "SELECT full_file_path from files where lfn='" + lname + "' and full_host_name = 'lustre'";

  odbc::Statement *stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet *rs = stmt->executeQuery(sqlquery);

  if (rs->next())
  {
    pfn = rs->getString(1);
    bret = true;
  }
  delete rs;
  delete stmt;
  return bret;
}

bool FROG::MinIOSearch(const std::string &lname)
{
  bool bret = false;
  if (!GetConnection())
  {
    return bret;
  }
  std::string sqlquery = "SELECT full_file_path from files where lfn='" + lname + "' and full_host_name = 'lustre'";

  if (Verbosity() > 1)
  {
    std::cout << "sql query:" << std::endl
	      << sqlquery << std::endl;
  }
  odbc::Statement *stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet *rs = stmt->executeQuery(sqlquery);

  if (rs->next())
  {
    pfn = rs->getString(1);
    std::string toreplace("/sphenix/lustre01/sphnxpro");
    size_t strpos = pfn.find(toreplace);
    if (strpos == std::string::npos)
    {
      std::cout << " could not locate " << toreplace 
		<< " in full file path " << pfn << std::endl;
      exit(1);
    }
    else if (strpos > 0)
    {
      std::cout << "full file path " << pfn 
		<< "does not start with " << toreplace << std::endl;
      exit(1);

    }
    pfn.replace(pfn.begin(),pfn.begin()+toreplace.size(),"s3://sphenixs3.rcf.bnl.gov:9000");
    bret = true;
  }
  delete rs;
  delete stmt;
  return bret;
}

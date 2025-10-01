#include "DBInterface.h"


#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <cassert>
#include <random> // For retrying connections
#include <sstream>
#include <string>
#include <chrono>
#include <iostream>
#include <thread>

#include <ffaobjects/RunHeader.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

DBInterface *DBInterface::__instance{nullptr};

DBInterface *DBInterface::instance()
{
  if (! __instance)
  {
    __instance = new DBInterface();
  }
  return __instance;
}
  
DBInterface::DBInterface(const std::string& name)
  : SubsysReco(name)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->addNewSubsystem(this);
}

int DBInterface::InitRun(PHCompositeNode* topNode) 
{
  if (Verbosity() > 1) { 
    std::cout << "Get run header" << std::endl;
  }
  
  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader) {
    std::cout << "can't find runheader" << std::endl;
    return 1;
  }
  std::cout << "run number: " << runheader->get_RunNumber() << std::endl;


  return 0;
}

int DBInterface::getRunTime(int runnumber) {
  odbc::Connection *m_OdbcConnection = getDBConnection("daq");

  std::string sql = "SELECT brtimestamp FROM run WHERE runnumber = " + std::to_string(runnumber) + ";";
  
  if (Verbosity() > 1) {
    std::cout << sql << std::endl;
  }
  odbc::Statement* stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet* resultSet = stmt->executeQuery(sql);

  if (resultSet && resultSet->next()) {
    odbc::Timestamp brtimestamp = resultSet->getTimestamp("brtimestamp");
    std::cout << brtimestamp.toString() << std::endl; // brtimestamp is in 'America/New_York' time zone
  }

  delete resultSet; 
  delete stmt; 
  delete m_OdbcConnection;
  return 0;

}

int DBInterface::End(PHCompositeNode* /*topNode*/)
{
  return 0;
}

double DBInterface::getDVal(const std::string &name)
{
  std::cout << name << std::endl;
  return 0.;
}

odbc::Connection *DBInterface::getDBConnection(const std::string &dbname, int verbosity)
{
  std::random_device ran_dev;
  std::seed_seq seeds {ran_dev(), ran_dev(), ran_dev()}; //...
  std::mt19937_64 mersenne_twister(seeds);
  std::uniform_int_distribution<> uniform(m_MIN_SLEEP_DUR, m_MAX_SLEEP_DUR);
  odbc::Connection* dbcon {nullptr};
  int total_slept_ms {0};
  int num_tries;
  for(num_tries = 0; num_tries < m_MAX_NUM_RETRIES; ++num_tries)
  {
    try
    {
      dbcon = odbc::DriverManager::getConnection(dbname, "", "");
    }
    catch (odbc::SQLException &e)
    {
      std::cerr << PHWHERE 
                << ": SQL Exception: "
                << e.getMessage() << std::endl;
    }

    if(dbcon)
    {
      ++num_tries;
      break;
    }
    int sleep_time_ms = uniform(mersenne_twister);
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time_ms));
    total_slept_ms += sleep_time_ms;

    if(0 < verbosity)
    {
      std::cout << PHWHERE
                << "Connection unsuccessful, Sleeping for addtional " << sleep_time_ms << " ms" << std::endl;
    }
  }
  if(1 < verbosity)
   {
     std::cout << PHWHERE 
               << ": Connection successful (" << num_tries << " attempts, " << total_slept_ms << " ms)" << std::endl;
   }
  if(!dbcon)
  {
    std::cerr << PHWHERE
              << ": DB Connection failed after " << m_MAX_NUM_RETRIES << " retries\n"
              << "Abandoning query" << std::endl;
    return nullptr;
  }
  return dbcon;
}



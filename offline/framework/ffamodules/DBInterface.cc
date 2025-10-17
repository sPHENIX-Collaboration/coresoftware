#include "DBInterface.h"

#include <ffaobjects/RunHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>  // For retrying connections
#include <sstream>
#include <string>
#include <thread>

DBInterface *DBInterface::__instance{nullptr};

DBInterface *DBInterface::instance()
{
  if (!__instance)
  {
    __instance = new DBInterface();
  }
  return __instance;
}

DBInterface::~DBInterface()
{
  if (!m_OdbcConnectionMap.empty())
  {
    for (auto iter : m_OdbcConnectionMap)
    {
      delete iter.second;
    }
    m_OdbcConnectionMap.clear();
  }
  if (!m_OdbcStatementMap.empty())
  {
    for (auto iter : m_OdbcStatementMap)
    {
      delete iter.second;
    }
    m_OdbcStatementMap.clear();
  }
  return;
}

DBInterface::DBInterface(const std::string &name)
  : SubsysReco(name)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->addNewSubsystem(this);
}

int DBInterface::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Get run header" << std::endl;
  }

  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader)
  {
    std::cout << "can't find runheader" << std::endl;
    return 1;
  }
  std::cout << "run number: " << runheader->get_RunNumber() << std::endl;

  return 0;
}

int DBInterface::process_event(PHCompositeNode * /*topNode*/)
{
  if (!m_OdbcConnectionMap.empty())
  {
    for (auto iter : m_OdbcConnectionMap)
    {
      delete iter.second;
    }
    m_OdbcConnectionMap.clear();
  }
  if (!m_OdbcStatementMap.empty())
  {
    // for (auto iter : m_OdbcStatementMap)
    // {
    //   delete iter.second;
    // }
    m_OdbcStatementMap.clear();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

double DBInterface::getDVal(const std::string &name)
{
  std::cout << name << std::endl;
  return 0.;
}

odbc::Connection *DBInterface::getDBConnection(const std::string &dbname, int verbosity)
{
  auto coniter = m_OdbcConnectionMap.find(dbname);
  if (coniter != m_OdbcConnectionMap.end())
  {
    return coniter->second;
  }
  std::random_device ran_dev;
  std::seed_seq seeds{ran_dev(), ran_dev(), ran_dev()};  //...
  std::mt19937_64 mersenne_twister(seeds);
  std::uniform_int_distribution<> uniform(m_MIN_SLEEP_DUR, m_MAX_SLEEP_DUR);
  odbc::Connection *dbcon{nullptr};
  int total_slept_ms{0};
  int num_tries;
  for (num_tries = 0; num_tries < m_MAX_NUM_RETRIES; ++num_tries)
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

    if (dbcon)
    {
      ++num_tries;
      break;
    }
    int sleep_time_ms = uniform(mersenne_twister);
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time_ms));
    total_slept_ms += sleep_time_ms;

    if (0 < verbosity)
    {
      std::cout << PHWHERE
                << "Connection unsuccessful, Sleeping for addtional " << sleep_time_ms << " ms" << std::endl;
    }
  }
  if (1 < verbosity)
  {
    std::cout << PHWHERE
              << ": Connection successful (" << num_tries << " attempts, " << total_slept_ms << " ms)" << std::endl;
  }
  if (!dbcon)
  {
    std::cout << PHWHERE
              << ": DB Connection failed after " << m_MAX_NUM_RETRIES << " retries\n"
              << "Abandoning query" << std::endl;
    return nullptr;
  }
  m_OdbcConnectionMap.insert(std::make_pair(dbname, dbcon));
  return dbcon;
}

odbc::Statement *DBInterface::getStatement(const std::string &dbname, int verbosity)
{
  auto statiter = m_OdbcStatementMap.find(dbname);
  if (statiter != m_OdbcStatementMap.end())
  {
    return statiter->second;
  }
  odbc::Connection *dbcon = getDBConnection(dbname, verbosity);
  odbc::Statement *statement = dbcon->createStatement();
  m_OdbcStatementMap.insert(std::make_pair(dbname, statement));
  return statement;
}

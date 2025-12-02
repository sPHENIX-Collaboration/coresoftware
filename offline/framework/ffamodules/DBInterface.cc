#include "DBInterface.h"

#include <sphenixodbc/ODBCInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/phool.h>

#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>  // For retrying connections
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
  delete m_ODBC;
  __instance = nullptr;
  return;
}

DBInterface::DBInterface(const std::string &name)
  : SubsysReco(name)
{
  m_ODBC = ODBCInterface::instance();
  Fun4AllServer *se = Fun4AllServer::instance();
  se->addNewSubsystem(this);
}

int DBInterface::process_event(PHCompositeNode * /*topNode*/)
{
  m_ODBC->Disconnect();
  return Fun4AllReturnCodes::EVENT_OK;
}

odbc::Connection *DBInterface::getDBConnection(const std::string &dbname)
{
  return m_ODBC->getDBConnection(dbname);
}

odbc::Statement *DBInterface::getStatement(const std::string &dbname)
{
  return m_ODBC->getStatement(dbname);
}

int DBInterface::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "Number of connection attempts" << std::endl;
    for (auto const &iter: m_ODBC->getNumConnection())
    {
      std::cout << "db: " << iter.first << ", attempts: " << iter.second << std::endl;
    }
    std::cout << "Total time slept: " << m_ODBC->SleepMS() << " ms" << std::endl;
    std::cout << "Total number of connection re-tries: " << m_ODBC->ConnectionTries() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void DBInterface::Print(const std::string & /*what*/) const
{
  if (m_ODBC->getNumConnection().empty())
  {
    std::cout << "No ODBC connections cached" << std::endl;
  }
  else
  {
    std::cout << "Odbc Connections Opened: " << std::endl;
    for (auto const &iter: m_ODBC->getNumConnection())
    {
      std::cout << "db: " << iter.first << ", opened: " << iter.second << std::endl;
    }
  }
  if (m_ODBC->getNumStatementUse().empty())
  {
    std::cout << "No odbc::Statements cached" << std::endl;
  }
  else
  {
    std::cout << "Odbc Statement use: " << std::endl;
    for (auto const &iter: m_ODBC->getNumStatementUse())
    {
      std::cout << "db: " << iter.first << ", Statement use: " << iter.second << std::endl;
    }
  }
  return;
}

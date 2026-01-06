#include "InttOdbcQuery.h"

#include <fun4all/DBInterface.h>

#include <phool/phool.h>

#include <odbc++/resultset.h>
#include <odbc++/statement.h>

#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>

int InttOdbcQuery::Query(int runnumber)
{
  m_query_successful = false;

// statement will be deleted by DBInterface
  odbc::Statement* statement = DBInterface::instance()->getStatement("daq");

  int iret = 0;
  iret += (QueryStreaming(statement, runnumber) != 0);
  //...

  m_query_successful = (iret == 0);

  return iret;
}

int InttOdbcQuery::QueryStreaming(odbc::Statement *statement, int runnumber)
{
  std::unique_ptr<odbc::ResultSet> result_set;

  std::string sched_data;
  try
  {
    std::string sql = "SELECT sched_data FROM gtm_scheduler WHERE vgtm=1 AND sched_entry = 1 AND runnumber = " + std::to_string(runnumber) + ";";
    result_set = std::unique_ptr<odbc::ResultSet>(statement->executeQuery(sql));
    if (result_set && result_set->next())
    {
      sched_data = result_set->getString("sched_data");
    }
  }
  catch (odbc::SQLException& e)
  {
    std::cerr << PHWHERE << "\n"
              << "\tSQL Exception:\n"
              << "\t" << e.getMessage() << std::endl;
    return 1;
  }

  if (std::string{"{17,55,24,54}"} == sched_data)
  {
    /// Streaming
    m_is_streaming = true;
  }
  else if (std::string{"{0,54,91,53}"} == sched_data)
  {
    /// Triggered
    m_is_streaming = false;
  }
  else
  {
    std::cerr << PHWHERE << "\n"
              << "\tUnexpected value for sched_data: '" << sched_data << "'" << std::endl;
    return 1;
  }

  if (m_verbosity)
  {
    std::cout << "\tsched_data: '" << sched_data << "'" << std::endl;
  }

  return 0;
}

int InttOdbcQuery::QueryType(odbc::Statement *statement, int runnumber)
{
  std::unique_ptr<odbc::ResultSet> result_set;
  m_type = "";

  try
  {
    std::string sql = "SELECT runtype FROM run WHERE runnumber = " + std::to_string(runnumber) + ";";
    result_set = std::unique_ptr<odbc::ResultSet>(statement->executeQuery(sql));
    if (result_set && result_set->next())
    {
      m_type = result_set->getString("runtype");
    }
  }
  catch (odbc::SQLException& e)
  {
    std::cerr << PHWHERE << "\n"
              << "\tSQL Exception:\n"
              << "\t" << e.getMessage() << std::endl;
    return 1;
  }

  if (m_verbosity)
  {
    std::cout << "\trun type: " << m_type << std::endl;
  }

  return 0;
}

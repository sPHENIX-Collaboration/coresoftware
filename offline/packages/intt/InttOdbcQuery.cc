/// Originally imitated from coresoftware/offline/database/rundb/CaloTime.*

#include "InttOdbcQuery.h"

#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>

#include <iostream>
#include <sstream>

#include <filesystem>
#include <random> // For retrying connections
#include <chrono> // For retrying connections
#include <thread> // For retrying connections

int InttOdbcQuery::Query(int runnumber)
{
  m_query_successful = false;

  std::random_device ran_dev;
  std::seed_seq seeds {ran_dev(), ran_dev(), ran_dev()}; //...
  std::mt19937_64 mersenne_twister(seeds);
  std::uniform_int_distribution<> uniform(m_MIN_SLEEP_DUR, m_MAX_SLEEP_DUR);

  int num_tries = 0;
  int total_slept_ms = 0;

  odbc::Connection* dbcon = nullptr;
  for(num_tries = 0; num_tries < m_MAX_NUM_RETRIES; ++num_tries)
  {
    try
    {
      dbcon = odbc::DriverManager::getConnection("daq", "", "");
    }
    catch (odbc::SQLException &e)
    {
      std::cerr << PHWHERE << "\n"
                << "\tSQL Exception:\n"
                << "\t" << e.getMessage() << std::endl;
    }

    if(dbcon)
    {
      ++num_tries;
      break;
    }

    int sleep_time_ms = uniform(mersenne_twister);
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time_ms));
    total_slept_ms += sleep_time_ms;

    if(1 < m_verbosity)
    {
      std::cout << PHWHERE << "\n"
                << "\tConnection unsuccessful\n"
                << "\tSleeping for addtional " << sleep_time_ms << " ms" << std::endl;
    }
  }

  if(m_verbosity)
  {
    std::cout << PHWHERE << "\n"
              << "\tConnection successful (" << num_tries << " attempts, " << total_slept_ms << " ms)" << std::endl;
  }

  if(!dbcon)
  {
    std::cerr << PHWHERE << "\n"
              << "\tDB Conntion failed after " << m_MAX_NUM_RETRIES << " retries\n"
              << "\tAbandoning query" << std::endl;
    return 1;
  }

  odbc::Statement* statement = dbcon->createStatement();

  int iret = 0;
  iret += (QueryStreaming(statement, runnumber) != 0);
  //...

  m_query_successful = (iret == 0);

  delete statement;
  delete dbcon;

  return iret;
}

int InttOdbcQuery::QueryStreaming(void* statement, int runnumber)
{
  odbc::ResultSet* result_set = nullptr;

  std::string sched_data;
  try
  {
    std::string sql = "SELECT sched_data FROM gtm_scheduler WHERE vgtm=1 AND sched_entry = 1 AND runnumber = " + std::to_string(runnumber) + ";";
    result_set = ((odbc::Statement*)statement)->executeQuery(sql);
    if(result_set && result_set->next())
    {
      sched_data = result_set->getString("sched_data");
    }
  }
  catch (odbc::SQLException &e)
  {
    std::cerr << PHWHERE << "\n"
              << "\tSQL Exception:\n"
              << "\t" << e.getMessage() << std::endl;
    delete result_set;
    return 1;
  }
  delete result_set;

  if(std::string{"{17,55,24,54}"} == sched_data)
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

  if(m_verbosity)
  {
    std::cout << "\tsched_data: '" << sched_data << "'" << std::endl;
  }

  return 0;
}

int InttOdbcQuery::QueryType(void* statement, int runnumber)
{
  odbc::ResultSet* result_set = nullptr;
  m_type = "";

  try
  {
    std::string sql = "SELECT runtype FROM run WHERE runnumber = " + std::to_string(runnumber) + ";";
    result_set = ((odbc::Statement*)statement)->executeQuery(sql);
    if(result_set && result_set->next())
    {
      m_type = result_set->getString("runtype");
    }
  }
  catch (odbc::SQLException &e)
  {
    std::cerr << PHWHERE << "\n"
              << "\tSQL Exception:\n"
              << "\t" << e.getMessage() << std::endl;
    delete result_set;
    return 1;
  }
  delete result_set;

  if(m_verbosity)
  {
    std::cout << "\trun type: " << m_type << std::endl;
  }

  return 0;
}

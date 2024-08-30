/// Originally imitated from coresoftware/offline/database/rundb/CaloTime.*

#include "InttOdbcQuery.h"

#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>

#include <iostream>
#include <sstream>

/// For retrying connections
#include <random>
#include <chrono>
#include <thread>

int InttOdbcQuery::Query(int runnumber)
{
  m_query_successful = false;

  /// Replace the following block with Singleton accessor? Something like:
  /// odbc::Connection* dbcon = ODBInterface::instance()->getConnection();
  /// Then ODBInterface implements the following 
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
  /// Replace the above block with Singleton accessor? Something like:
  /// odbc::Connection* dbcon = ODBInterface::instance()->getConnection();
  /// Then ODBInterface implements the above

  if(!dbcon)
  {
    std::cerr << PHWHERE << "\n"
              << "\tDB Conntion failed after " << m_MAX_NUM_RETRIES << " retries\n"
              << "\tAbandoning query" << std::endl;
    return 1;
  }

  std::string sql = "select sched_data from  gtm_scheduler where vgtm=1 and sched_entry = 1 and runnumber = " + std::to_string(runnumber) + ";";

  odbc::Statement *statement = dbcon->createStatement();
  odbc::ResultSet *result_set = statement->executeQuery(sql);
  std::string sched_data;

  if(result_set && result_set->next())
  {
	try
    {
      /// This can throw an odbc::SQLException if, for example, the column name is mistyped
      /// (you'll never guess how I figured that one out...)
      sched_data = result_set->getString("sched_data");
    }
    catch (odbc::SQLException &e)
    {
      std::cerr << PHWHERE << "\n"
                << "\tSQL Exception:\n"
                << "\t" << e.getMessage() << std::endl;
      delete result_set;
      delete statement;
      delete dbcon;
    
      return 1;
    }
  }

  if(m_verbosity)
  {
    std::cout << PHWHERE << "\n"
              << "\tsched_data: '" << sched_data << "'" << std::endl;
  }

  if(std::string{"{17,55,24,54}"} == sched_data)
  {
    /// Streaming
    m_query_successful = true;
    m_is_streaming = true;
  }
  else if (std::string{"{0,54,91,53}"} == sched_data)
  {
    /// Triggered
    m_query_successful = true;
    m_is_streaming = false;
  }

  delete result_set;
  delete statement;
  delete dbcon;

  return 0;
}

bool InttOdbcQuery::IsStreaming()
{
  if(!m_query_successful)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCall is not preceeded by successful database connection and query\n"
              << "\tValue returned will not necessarily be correct" << std::endl;
  }

  return m_is_streaming;
}

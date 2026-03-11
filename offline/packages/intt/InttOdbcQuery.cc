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

  m_intt_dac_values.clear();

  int iret = 0;
  iret += (QueryStreaming(statement, runnumber) != 0);
  iret += (QueryAllDACValues(statement, runnumber) != 0);
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

int InttOdbcQuery::QuerySingleDACValue(odbc::Statement *statement, int runnumber, int adc_value, int &DAC_value) // adc_value should be from 0 to 7
{
  std::unique_ptr<odbc::ResultSet> result_set;
  std::string column_name = "dac" + std::to_string(adc_value);
  DAC_value = -1;
  

  try
  {
    std::string sql = "SELECT " + column_name + " From intt_setting WHERE runnumber = " + std::to_string(runnumber) + ";";
    result_set = std::unique_ptr<odbc::ResultSet>(statement->executeQuery(sql));
    if (result_set && result_set->next())
    {
      DAC_value = result_set->getInt(column_name.c_str());
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
    std::cout << column_name << " of run " << runnumber <<" is " << DAC_value << std::endl;
  }

  return 0;
}

int InttOdbcQuery::QueryAllDACValues(odbc::Statement *statement, int runnumber)
{
  int error_count = 0;

  for (int i = 0; i < 8; ++i)
  {
    int DAC_value = -1;

    error_count += QuerySingleDACValue(statement, runnumber, i, DAC_value);

    m_intt_dac_values.push_back(DAC_value);
  }

  if (error_count != 0 || m_intt_dac_values.size() != 8 || std::find(m_intt_dac_values.begin(), m_intt_dac_values.end(), -1) != m_intt_dac_values.end())
  {
    std::cerr << PHWHERE << "\n"
              << "\tError retrieving DAC values. error_count: " << error_count
              << ", size of m_intt_dac_values: " << m_intt_dac_values.size() << std::endl;
    std:: cout << "In the dac map: ";
    for (size_t i = 0; i < m_intt_dac_values.size(); ++i)
    {
      std::cout << "dac" << i << ": " << m_intt_dac_values[i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "We then use the default mapping of intt_dac_values: {30, 45, 60, 90, 120, 150, 180, 210} for this run, runnumber : "<< runnumber << std::endl;
    m_intt_dac_values = default_intt_dac_values;

    return error_count;
  }

  return 0;
}
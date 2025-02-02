#include "MvtxRawDefs.h"

#include <odbc++/connection.h> // odbc::Connection
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <odbc++/drivermanager.h>
#pragma GCC diagnostic pop
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

#include <TRandom3.h>

#include <limits>
#include <string>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <thread>

uint8_t MvtxRawDefs::getStaveIndex( const uint8_t& lyrId, const uint8_t& stvId )
{
  return firstStaveIndex[lyrId] + stvId;
};

std::pair<uint8_t, uint8_t> const& MvtxRawDefs::get_flx_endpoint( const uint8_t& lyrId, const uint8_t& stvId )
{
  return stave_felix_map.at( getStaveIndex(lyrId, stvId) );
};

MvtxRawDefs::linkId_t MvtxRawDefs::decode_feeid( const uint16_t feeid )
{
  linkId_t ret = {};
// the static_cast< uint16_t> is needed to because the result of (feeid >> 12U)
// is promoted to int which then triggers a (correct) clang-tidy warning that
// a bitwise operation is performed on a signed integer
  ret.layer = static_cast<uint16_t>(feeid >> 12U) & 0x7U;
  ret.stave = feeid & 0x1FU;
  ret.gbtid = static_cast<uint16_t>(feeid >> 8U) & 0x3U;
  return ret;
};

// anonymous namespace
namespace {
  float getStrobeLengthFromOCDB(const int& runNumber)
  {
    constexpr uint8_t MAX_NUM_RETRIES = 10;
    std::unique_ptr<odbc::Connection> m_OdbcConnection(nullptr);
    int num_tries = 0;
    TRandom3 rnd = TRandom3();
    do
    {
      try
      {
        std::unique_ptr<odbc::Connection> temp(odbc::DriverManager::getConnection("mvtx_read","",""));
        m_OdbcConnection = std::move(temp);
      }
      catch (odbc::SQLException &e)
      {
        std::cout << " Exception caught during DriverManager::getConnection" << std::endl;
        std::cout << "Message: " << e.getMessage() << std::endl;
      }
      ++num_tries;
      int wait = (int) 20 + rnd.Uniform()*300;
      if (!m_OdbcConnection)
      {
        std::this_thread::sleep_for(std::chrono::seconds(wait));  // sleep 30 seconds before retry
      }
    }
    while ((! m_OdbcConnection) && (num_tries < MAX_NUM_RETRIES));

    if (! m_OdbcConnection)
    {
      return std::numeric_limits<float>::quiet_NaN();
    }

    std::string sqlcmd = "SELECT strobe FROM mvtx_strobe_offline WHERE runnumber = " + std::to_string(runNumber) + ";";

    std::unique_ptr<odbc::Statement> statement (m_OdbcConnection->createStatement());
    std::unique_ptr<odbc::ResultSet> resultSet (statement->executeQuery(sqlcmd));

    if (resultSet && resultSet->next())
    {
      float strobe = resultSet->getFloat("strobe");
      std::cout << "MVTX strobe read from ocdb: " << strobe << std::endl;
      return strobe;
    }

    return std::numeric_limits<float>::quiet_NaN();
  }

  float getStrobeLengthFromDAQ(const int& runNumber)
  {
    float strobeWidth = std::numeric_limits<float>::quiet_NaN();

    std::string executable_command = "psql -h sphnxdaqdbreplica daq --csv -c \"SELECT strobe FROM mvtx_strobe WHERE hostname = \'mvtx0\' AND runnumber = ";
    executable_command += std::to_string(runNumber);
    executable_command += ";\" | tail -n 1";

    std::array<char, 128> buffer = {};
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(executable_command.c_str(), "r"), pclose);
    if (!pipe)
    {
      throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
    {
      result += buffer.data();
    }
    try
    {
      strobeWidth = stof(result);
    }
    catch (std::invalid_argument const& ex)
    {
      std::cout << "mvtx_utils::getStrobeLength() Run number " << runNumber << " has no strobe length in the DAQ database, returning NAN" << std::endl;
    }

    std::cout << "MVTX strobe read from daq db: " << strobeWidth << std::endl;
    return strobeWidth;
  }
}

float MvtxRawDefs::getStrobeLength(const int& runNumber)
{
  {
    float strobe_length = getStrobeLengthFromOCDB(runNumber);
    return (strobe_length != std::numeric_limits<float>::quiet_NaN()) ? \
            strobe_length : getStrobeLengthFromDAQ(runNumber);
  }
}
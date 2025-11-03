#include "MvtxRawDefs.h"

#include <ffamodules/DBInterface.h>

#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

uint8_t MvtxRawDefs::getStaveIndex(const uint8_t& lyrId, const uint8_t& stvId)
{
  return firstStaveIndex[lyrId] + stvId;
};

std::pair<uint8_t, uint8_t> const& MvtxRawDefs::get_flx_endpoint(const uint8_t& lyrId, const uint8_t& stvId)
{
  return stave_felix_map.at(getStaveIndex(lyrId, stvId));
};

MvtxRawDefs::linkId_t MvtxRawDefs::decode_feeid(const uint16_t feeid)
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
namespace
{
  float getStrobeLengthFromOCDB(const int& runNumber)
  {
    odbc::Statement *statement = DBInterface::instance()->getStatement("mvtx_read");

    std::string sqlcmd = "SELECT strobe FROM mvtx_strobe_offline WHERE runnumber = " + std::to_string(runNumber);

    std::unique_ptr<odbc::ResultSet> resultSet(statement->executeQuery(sqlcmd));

    if (resultSet && resultSet->next())
    {
      float strobe = resultSet->getFloat("strobe");
      return strobe;
    }

    return std::numeric_limits<float>::quiet_NaN();
  }

  float getStrobeLengthFromDAQ(const int& runNumber)
  {
    float strobeWidth = std::numeric_limits<float>::quiet_NaN();
    odbc::Statement *statement = DBInterface::instance()->getStatement("daq");
    std::string sqlcmd = "SELECT strobe FROM mvtx_strobe WHERE hostname = 'mvtx0' AND runnumber = " + std::to_string(runNumber);
    std::unique_ptr<odbc::ResultSet> resultSet(statement->executeQuery(sqlcmd));
    if (resultSet && resultSet->next())
    {
      strobeWidth = resultSet->getFloat("strobe");
    }
    return strobeWidth;
  }
}  // namespace

float MvtxRawDefs::getStrobeLength(const int& runNumber)
{
  {
    float strobe_length = getStrobeLengthFromOCDB(runNumber);
    return (!std::isnan(strobe_length)) ? strobe_length : getStrobeLengthFromDAQ(runNumber);
  }
}

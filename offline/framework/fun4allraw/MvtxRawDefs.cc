#include "MvtxRawDefs.h"

#include <limits>
#include <string>
#include <memory>
#include <stdexcept>
#include <iostream>
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



float MvtxRawDefs::getStrobeLength(const int& runNumber)
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
  
  return strobeWidth;
}
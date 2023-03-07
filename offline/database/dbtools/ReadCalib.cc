#include "ReadCalib.h"

#include <phool/recoConsts.h>

#include <nlohmann/json.hpp>  // every method returns a json object
#include <sphenixnpc.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

int ReadCalib::listGlobalTags()
{
  std::cout << sphenixnpc::getGlobalTags() << std::endl;
  return 0;
}

int ReadCalib::listPayloadTypes()
{
  std::cout << sphenixnpc::getPayloadTypes() << std::endl;
  return 0;
}

std::string ReadCalib::getCalibrationFile(const std::string &type, uint64_t iov) const
{
  recoConsts *rc = recoConsts::instance();
  std::string globaltag;
  if (rc->FlagExist("CDB_GLOBALTAG"))
  {
    globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  }
  else
  {
    std::cout << "rc flag for globat tag CDB_GLOBALTAG not set" << std::endl;
    exit(1);
  }
  return getCalibrationFile(globaltag, type, iov);
}

std::string ReadCalib::getCalibrationFile(const std::string &globaltag, const std::string &type, uint64_t iov) const
{
  nlohmann::json result = sphenixnpc::get(globaltag, type, iov);
  return result.at("msg");
}

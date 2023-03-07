#include "InsertCalib.h"

#include <phool/recoConsts.h>

#include <sphenixnpc.hpp>

#include <nlohmann/json.hpp>  // every method returns a json object

#include <iostream>
#include <stdexcept>  // for out_of_range
#include <string>

int InsertCalib::createGlobalTag(const std::string &tag)
{
  std::cout << sphenixnpc::createGlobalTag(tag) << std::endl;
  return 0;
}
int InsertCalib::listGlobalTags()
{
  std::cout << sphenixnpc::getGlobalTags() << std::endl;
  nlohmann::json result = sphenixnpc::getGlobalTags();

  std::cout << "code:" << result.at("code") << std::endl;
  return 0;
}
int InsertCalib::deleteGlobalTag(const std::string &tag)
{
  std::cout << sphenixnpc::deleteGlobalTag(tag) << std::endl;
  return 0;
}

int InsertCalib::createPayloadType(const std::string &payloadtype)
{
  std::cout << sphenixnpc::createPayloadType(payloadtype) << std::endl;
  return 0;
}

int InsertCalib::listPayloadTypes()
{
  std::cout << sphenixnpc::getPayloadTypes() << std::endl;
  return 0;
}

int InsertCalib::insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start)
{
  recoConsts *rc = recoConsts::instance();
  std::string globaltag;
  if (rc->FlagExist("CDB_GLOBALTAG"))
  {
    globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  }
  std::cout << sphenixnpc::createGlobalTag(globaltag) << std::endl;
  std::cout << sphenixnpc::createPayloadType(payloadtype) << std::endl;
  ;
  std::cout << "inserting " << fname << std::endl;
  std::cout << sphenixnpc::insertPayload(globaltag, payloadtype, fname, iov_start) << std::endl;
  return 0;
}

int InsertCalib::insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start, uint64_t iov_end)
{
  recoConsts *rc = recoConsts::instance();
  std::string globaltag;
  if (rc->FlagExist("CDB_GLOBALTAG"))
  {
    globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  }
  std::cout << sphenixnpc::createGlobalTag(globaltag) << std::endl;
  std::cout << sphenixnpc::createPayloadType(payloadtype) << std::endl;
  ;
  std::cout << "inserting " << fname << std::endl;
  std::cout << sphenixnpc::insertPayload(globaltag, payloadtype, fname, iov_start, iov_end) << std::endl;
  return 0;
}

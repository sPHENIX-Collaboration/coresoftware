#include "CDBUtils.h"

#include "SphenixClient.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <map>
#include <set>
#include <stdexcept>  // for out_of_range
#include <tuple>

CDBUtils::CDBUtils()
  : cdbclient(new SphenixClient())
{
}

CDBUtils::CDBUtils(const std::string &globaltag)
  : cdbclient(new SphenixClient(globaltag))
{
}

int CDBUtils::createGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->createGlobalTag(tagname);
  int iret = resp["code"];
  nlohmann::json msgcont = resp["msg"];
  std::cout << msgcont << std::endl;
  return iret;
}

int CDBUtils::deleteGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->deleteGlobalTag(tagname);
  int iret = resp["code"];
  nlohmann::json msgcont = resp["msg"];
  std::cout << msgcont << std::endl;
  return iret;
}

int CDBUtils::lockGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->lockGlobalTag(tagname);
  int iret = resp["code"];
  nlohmann::json msgcont = resp["msg"];
  std::cout << msgcont << std::endl;
  return iret;
}

int CDBUtils::unlockGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->unlockGlobalTag(tagname);
  int iret = resp["code"];
  nlohmann::json msgcont = resp["msg"];
  std::cout << "message: " << msgcont << std::endl;
  return iret;
}

void CDBUtils::clearCache()
{
  nlohmann::json resp = cdbclient->clearCache();
  std::cout << resp["msg"] << std::endl;
  return;
}

std::string CDBUtils::getUrl(const std::string &type, uint64_t iov)
{
  nlohmann::json resp = cdbclient->getUrl(type, iov);
  return resp["msg"];
}

int CDBUtils::createPayloadType(const std::string &pt)
{
  return cdbclient->createDomain(pt);
}

void CDBUtils::listPayloadIOVs(uint64_t iov)
{
  nlohmann::json resp = cdbclient->getPayloadIOVs(iov);
  if (resp["code"] != 0)
  {
    std::cout << resp["msg"] << std::endl;
    return;
  }
  nlohmann::json payload_iovs = resp["msg"];
  std::map<std::string, std::tuple<std::string, uint64_t, uint64_t>> iovs;
  for (auto &[pt, val] : payload_iovs.items())
  {
    std::string url = val["payload_url"];
    uint64_t bts = val["minor_iov_start"];
    uint64_t ets = val["minor_iov_end"];
    iovs.insert(std::make_pair(pt, std::make_tuple(url, bts, ets)));
  }
  for (const auto &it : iovs)
  {
    std::cout << it.first << ": " << std::get<0>(it.second)
              << ", begin ts: " << std::get<1>(it.second)
              << ", end ts: " << std::get<2>(it.second)
              << std::endl;
  }
  return;
}

int CDBUtils::cloneGlobalTag(const std::string &source, const std::string &target)
{
  nlohmann::json resp = cdbclient->getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  std::set<std::string> gtset;
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    gtset.insert(exist_gt);
  }
  if (gtset.find(source) == gtset.end())
  {
    std::cout << "source tag " << source << " does not exist" << std::endl;
    return -1;
  }
  if (gtset.find(target) != gtset.end())
  {
    std::cout << "Target tag " << target << " exists, delete it first" << std::endl;
    return -1;
  }
  resp = cdbclient->cloneGlobalTag(source, target);
  int iret = resp["code"];
  std::cout << resp["msg"] << std::endl;
  return iret;
}

void CDBUtils::listGlobalTags()
{
  nlohmann::json resp = cdbclient->getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  std::set<std::string> globaltags;
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    globaltags.insert(exist_gt);
  }
  for (const auto &it : globaltags)
  {
    std::cout << "global tag: " << it << std::endl;
  }
  return;
}

void CDBUtils::listPayloadTypes()
{
  nlohmann::json resp = cdbclient->getPayloadTypes();
  nlohmann::json msgcont = resp["msg"];
  std::set<std::string> payloadtypes;
  for (auto &it : msgcont.items())
  {
    std::string exist_pl = it.value().at("name");
    payloadtypes.insert(exist_pl);
  }
  for (const auto &it : payloadtypes)
  {
    std::cout << "payload type: " << it << std::endl;
  }
  return;
}

int CDBUtils::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start)
{
  if (!isGlobalTagSet())
  {
    std::cout << "No Global Tag set" << std::endl;
    return -1;
  }
  nlohmann::json resp = cdbclient->insertPayload(pl_type, file_url, iov_start);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error inserting payload " << file_url << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp["msg"] << std::endl;
  }
  return iret;
}

int CDBUtils::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end)
{
  if (!isGlobalTagSet())
  {
    std::cout << "No Global Tag set" << std::endl;
    return -1;
  }
  nlohmann::json resp = cdbclient->insertPayload(pl_type, file_url, iov_start, iov_end);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error inserting payload " << file_url << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp["msg"] << std::endl;
  }
  return iret;
}

int CDBUtils::setGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->setGlobalTag(tagname);
  int iret = resp["code"];
  std::cout << "message: " << resp["msg"] << std::endl;
  return iret;
}

bool CDBUtils::isGlobalTagSet()
{
  return cdbclient->isGlobalTagSet();
}

void CDBUtils::Verbosity(int i)
{
  if (cdbclient)
  {
    cdbclient->Verbosity(i);
  }
  m_Verbosity = i;
}

int CDBUtils::deletePayloadIOV(const std::string &pl_type, uint64_t iov_start)
{
  nlohmann::json resp = cdbclient->deletePayloadIOV(pl_type, iov_start);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error deleting payload iov, type " << pl_type
              << ", iov_start: " << iov_start
              << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp["msg"] << std::endl;
  }
  return iret;
}

int CDBUtils::deletePayloadIOV(const std::string &pl_type, uint64_t iov_start, uint64_t iov_end)
{
  nlohmann::json resp = cdbclient->deletePayloadIOV(pl_type, iov_start, iov_end);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error deleting payload iov, type " << pl_type
              << ", iov_start: " << iov_start
              << ", iov_end: " << iov_end
              << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp["msg"] << std::endl;
  }
  return iret;
}
